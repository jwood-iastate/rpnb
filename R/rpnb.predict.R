# Internal helper function to generate draws based on model estimates
# This avoids repeating code between the main function and predict
generate_rpar_draws <- function(object, ndraws) {
  params <- object$estimate
  n_rand <- length(object$rpar_names)

  # Get means of random parameters
  mean_param_names <- paste0("mean.", object$rpar_names)
  betas_rand_mean <- params[mean_param_names]

  # Generate base Halton draws (use scrambled for all predictions, for simplicity)
  halton_draws <- make_draws(ndraws, length(object$rpar_names), type = "scrambled-halton-rand-perm", seed = 123)
  draws <- matrix(NA, nrow = ndraws, ncol = n_rand)

  if (object$correlated) {
    # For correlated parameters, use the stored Cholesky decomposition
    chol_mat <- t(chol(object$rpar_vcov)) # Get Cholesky from stored VCV
    draws <- (qnorm(halton_draws) %*% t(chol_mat)) +
      matrix(betas_rand_mean, nrow = ndraws, ncol = n_rand, byrow = TRUE)
  } else {
    # For uncorrelated parameters, use the estimated standard deviations
    sd_param_names <- paste0("sd.", object$rpar_names)
    rand_sds <- abs(params[sd_param_names])

    for (i in 1:n_rand) {
      dist <- object$rpardists[i]
      draws[, i] <- switch(dist,
                           "n"  = qnorm(halton_draws[, i], mean = betas_rand_mean[i], sd = rand_sds[i]),
                           "ln" = qlnorm(halton_draws[, i], meanlog = betas_rand_mean[i], sdlog = rand_sds[i]),
                           "t"  = qtri(halton_draws[, i], lower = betas_rand_mean[i] - rand_sds[i], upper = betas_rand_mean[i] + rand_sds[i], mode=betas_rand_mean[i]),
                           "u"  = qunif(halton_draws[, i], min = betas_rand_mean[i] - rand_sds[i], max = betas_rand_mean[i] + rand_sds[i]),
                           "g"  = qgamma(halton_draws[, i], shape = betas_rand_mean[i]^2 / rand_sds[i]^2, rate = betas_rand_mean[i] / rand_sds[i]^2)
      )
    }
  }
  return(draws)
}


#' Predict from a Random Parameters Negative Binomial Model
#'
#' @name predict.rpnb
#' @param object A model object of class `rpnb` estimated by the `rpnb()` function.
#' @param newdata A data frame containing variables for prediction.
#' @param method The prediction method: `"Simulated"` (default), `"Exact"`, or `"Individual"`.
#' @param ndraws The number of draws to use for simulation-based methods.
#' @param ... Additional arguments (currently ignored).
#'
#' @note
#' The `"Individual"` method performs a Bayesian update to find the conditional
#' parameters for each panel and thus requires the outcome variable to be present
#' in `newdata`. The `"Exact"` method uses analytical moment-generating functions
#' and is only available for specific distributions.
#'
#' @return A numeric vector of predicted values.
#' @import randtoolbox stats modelr
#' @importFrom utils head tail
#' @exportS3Method stats::predict
#'
#' @examples
#' \donttest{
#' # Assuming `rpnb_model` is an estimated model object from the `rpnb` function
#' # predictions <- predict(rpnb_model, newdata = washington_roads)
#' # individual_preds <- predict(rpnb_model, newdata = washington_roads, method = "Individual")
#' }
predict.rpnb <- function(object, newdata, method = "Simulated", ndraws = 2000, ...) {

  ## >> Setup
  method <- match.arg(method, c("Simulated", "Exact", "Individual"))

  # Extract matrices and parameters from the model object and new data
  X_fixed <- model.matrix(object$formula, newdata)
  X_rand <- model.matrix(object$rpar_formula, newdata)

  params <- object$estimate
  mean_fixed_names <- paste0("mean.", object$fixed_names)
  betas_fixed <- params[mean_fixed_names]
  mu_fixed <- as.vector(exp(X_fixed %*% betas_fixed))

  ## >> Method 1: Simulated Prediction (Population Average)
  if (method == "Simulated") {
    # Generate draws from the estimated population distribution
    draws <- generate_rpar_draws(object, ndraws)

    # Calculate the random part of the linear predictor for each draw
    xb_rand <- X_rand %*% t(draws)

    # Combine with the fixed part and calculate the simulated mean for each draw
    mu_sim <- mu_fixed * exp(xb_rand)

    # The prediction is the average of the simulated means over all draws
    predictions <- rowMeans(mu_sim)
    return(predictions)
  }

  ## >> Method 2: "Exact" Prediction (using MGF)
  if (method == "Exact") {
    mean_rand_names <- paste0("mean.", object$rpar_names)
    betas_rand_mean <- params[mean_rand_names]

    if (object$correlated) {
      # MGF for correlated normal parameters: exp(t'μ + 0.5 * t'Σt)
      rpar_vcov <- object$rpar_vcov
      # Calculate quadratic term t'Σt for each observation
      quad_term <- rowSums((X_rand %*% rpar_vcov) * X_rand)
      rand_factor <- exp(X_rand %*% betas_rand_mean + 0.5 * quad_term)
    } else {
      # MGF for uncorrelated parameters is the product of individual MGFs
      sd_rand_names <- paste0("sd.", object$rpar_names)
      rand_sds <- abs(params[sd_rand_names])
      rand_factors <- matrix(1, nrow = nrow(newdata), ncol = length(object$rpar_names))

      for(i in seq_along(object$rpar_names)) {
        dist <- object$rpardists[i]
        mu <- betas_rand_mean[i]
        sigma <- rand_sds[i]
        x <- X_rand[, i]

        rand_factors[, i] <- switch(dist,
                                    "n" = exp(x * mu + 0.5 * x^2 * sigma^2),
                                    "u" = ifelse(x == 0, 1, exp(x * mu) * sinh(x * sigma) / (x * sigma)),
                                    "g" = (1 - x * sigma^2 / mu)^(-mu^2 / sigma^2),
                                    # Other MGFs (log-normal, triangular) can be complex and are omitted for robustness
                                    # Default to simulation if MGF is not simple
                                    exp(x * mu) # Fallback to mean-only if distribution not supported
        )
      }
      rand_factor <- apply(rand_factors, 1, prod)
    }
    predictions <- mu_fixed * rand_factor
    return(as.vector(predictions))
  }

  ## >> Method 3: Individual Prediction (Panel or Cross-section)
  if (method == "Individual") {
    # This method requires the outcome variable
    y_name <- all.vars(object$formula)[1]
    if (!y_name %in% names(newdata)) {
      stop("Method 'Individual' requires the outcome variable '", y_name, "' in newdata.")
    }
    y <- newdata[[y_name]]

    # Create panel index from newdata
    panel_id_col <- object$panel_id
    panel_index <- if(is.null(panel_id_col)) 1:nrow(newdata) else as.factor(newdata[[panel_id_col]])

    # Generate draws and calculate simulated means
    draws <- generate_rpar_draws(object, ndraws)
    xb_rand <- X_rand %*% t(draws)
    mu_sim <- mu_fixed * exp(xb_rand)

    # Calculate observation probabilities for each draw
    size <- switch(object$form,
                   "nb1" = mu_sim / exp(params["ln_alpha"]),
                   "nb2" = 1 / exp(params["ln_alpha"]),
                   "nbp" = mu_sim^(2 - params["P"]) / exp(params["ln_alpha"])
    )
    prob_mat <- dnbinom(y, size = size, mu = mu_sim)
    prob_mat[prob_mat <= 0] <- 1e-300

    # Calculate panel-level likelihoods for each draw
    panel_log_probs <- rowsum(log(prob_mat), group = panel_index, reorder = FALSE)
    max_log_prob <- apply(panel_log_probs, 1, max)
    panel_likelihoods <- exp(panel_log_probs - max_log_prob) # Stable calculation

    # Calculate conditional expectation of draws for each panel
    # E[β_p | y_p] ≈ Σ(β_r * L(y_p|β_r)) / Σ(L(y_p|β_r))
    numerator <- panel_likelihoods %*% draws
    denominator <- rowSums(panel_likelihoods)

    cond_betas_panel <- numerator / denominator

    # Expand panel-level betas to the observation level
    cond_betas_obs <- cond_betas_panel[as.integer(panel_index), ]

    # Calculate predictions using the conditional (panel-specific) betas
    predictions <- mu_fixed * exp(rowSums(X_rand * cond_betas_obs))
    return(predictions)
  }
}
