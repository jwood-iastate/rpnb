#' Estimate a Random Parameters Negative Binomial Model
#'
#' @description
#' Estimates a random parameters negative binomial (RPNB) model using maximum
#' simulated likelihood. Supports NB-1, NB-2, and NB-P specifications,
#' correlated or uncorrelated random parameters, various parameter distributions,
#' and panel data structures.
#'
#' @name rpnb
#' @param formula A standard R formula for the model's fixed parameters.
#' @param rpar_formula A one-sided formula specifying the random parameters (e.g., `~ var1 + var2`).
#' @param data A data frame containing the variables for the model.
#' @param panel ## >> A character string specifying the column name in `data` that
#'   identifies the panel groups (e.g., individuals). If `NULL` (the default),
#'   each observation is treated as a separate panel (cross-sectional data).
#' @param form The negative binomial form: `"nb2"`, `"nb1"`, or `"nbp"`.
#' @param rpardists An optional named character vector specifying the distribution
#'   for each random parameter. Options: `"n"` (normal), `"ln"` (log-normal),
#'   `"t"` (triangular), `"u"` (uniform), `"g"` (gamma). Defaults to normal for all.
#' @param ndraws The number of Halton draws for simulation.
#' @param scrambled Logical. If `TRUE`, uses scrambled Halton draws.
#' @param correlated Logical. If `TRUE`, estimates correlated random parameters
#'   (forces all to be normally distributed).
#' @param method The optimization algorithm to be used by `maxLik::maxLik`.
#' @param max.iters Maximum number of iterations for the optimizer.
#' @param start.vals An optional named vector of starting values.
#' @param print.level Verbosity level for the optimization process (0, 1, or 2).
#'
#' @return A list object of class `rpnb` containing the model results, including
#'   coefficients, variance-covariance matrix, log-likelihood, and other model details.
#'
#' @import maxLik stats modelr
#' @importFrom MASS glm.nb
#' @importFrom utils head tail
#' @include tri.R
#'
#' @export
#' @examples
#' \donttest{
#'
#' ## Random Parameters Negative Binomial model (NB-1)
#' data("washington_roads")
#' nb1.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = FALSE,
#'                rpardists = c(intercept="u", speed50="t"),
#'                form = 'nb1',
#'                method = "bfgs",
#'                print.level = 2)
#'
#' summary(nb1.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-2)
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = TRUE,
#'                form = 'nb2',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nb2.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-P)
#' nbp.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                rpardists = c(intercept = "u", speed50 = "n"),
#'                ndraws = 100,
#'                correlated = FALSE,
#'                form = 'nbp',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nbp.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-2) with panel
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = TRUE,
#'                form = 'nb2',
#'                panel = "ID",
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nb2.rp)}
rpnb <- function(formula, rpar_formula, data, panel = NULL, form = "nb2",
                 rpardists = NULL, ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, print.level = 0) {

  ## >> Improvement: Added panel argument to validation
  stopifnot(
    inherits(data, "data.frame"),
    inherits(formula, "formula"),
    inherits(rpar_formula, "formula"),
    is.null(panel) || (is.character(panel) && length(panel) == 1 && panel %in% names(data)),
    is.numeric(ndraws) && ndraws > 0,
    is.logical(scrambled),
    is.logical(correlated)
  )
  form <- match.arg(form, c("nb1", "nb2", "nbp"))

  # Set up model matrices and response variable
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X_fixed <- model.matrix(formula, data)
  X_rand <- model.matrix(rpar_formula, data)
  rpar_names <- colnames(X_rand)
  fixed_names <- colnames(X_fixed)

  ## >> Improvement: Efficient base R panel index creation
  if (is.null(panel)) {
    # Cross-sectional: each observation is its own panel
    panel_index <- 1:nrow(data)
  } else {
    # Panel data: create a factor from the specified panel column
    panel_index <- as.factor(data[[panel]])
  }
  n_panels <- length(unique(panel_index))


  # (Checks for intercept, correlation, and rpardists remain the same)
  if ("(Intercept)" %in% fixed_names && "(Intercept)" %in% rpar_names) {
    stop("The intercept cannot be both a fixed and random parameter.")
  }
  if (correlated && length(rpar_names) < 2) {
    message("Correlation requires at least two random parameters. Setting `correlated = FALSE`.")
    correlated <- FALSE
  }
  if (correlated && !is.null(rpardists)) {
    message("When `correlated = TRUE`, all random parameters are assumed to be normal. Ignoring `rpardists`.")
    rpardists <- NULL
  }

  if (!correlated && is.null(rpardists)) {
    rpardists <- rep("n", length(rpar_names))
    names(rpardists) <- rpar_names
  }
  if (!is.null(rpardists)) {
    names(rpardists) <- gsub("intercept|constant", "(Intercept)", names(rpardists), ignore.case = TRUE)
    if (!setequal(names(rpardists), rpar_names)) {
      stop("Names in `rpardists` must exactly match the random parameter names.")
    }
    rpardists <- rpardists[rpar_names]
  }

  # Generate Halton draws
  if (scrambled){
    halton_draws <- make_draws(ndraws, length(rpar_names), type = "scrambled-halton-rand-perm", seed = 123)
  }else{
    halton_draws <- make_draws(ndraws, length(rpar_names), type = "standard-halton")
  }

  # Smart starting values
  if (is.null(start.vals)) {
    full_formula <- update(formula, paste("~ . +", paste(rpar_names[rpar_names != "(Intercept)"], collapse = "+")))
    nb_model <- tryCatch(MASS::glm.nb(full_formula, data = data), error = function(e) glm(full_formula, data=data, family=poisson))
    start_betas <- coef(nb_model)
    all_vars <- c(fixed_names, rpar_names)
    start <- rep(0, length(all_vars)); names(start) <- all_vars
    common_coefs <- intersect(names(start), names(start_betas))
    start[common_coefs] <- start_betas[common_coefs]

    # Define parameter names based on correlation status
    if (correlated) {
      chol_mat <- diag(0.1, nrow = length(rpar_names))
      chol_vals <- chol_mat[lower.tri(chol_mat, diag = TRUE)]
      chol_names <- paste0("chol.", 1:length(chol_vals))
      start <- c(start, chol_vals)
      names(start) <- c(paste0("mean.", all_vars), chol_names)
    } else {
      sd_names <- paste0("sd.", rpar_names)
      start <- c(start, rep(0.1, length(rpar_names)))
      names(start) <- c(paste0("mean.", all_vars), sd_names)
    }

    # Add dispersion and P parameters
    theta <- if(!is.null(nb_model$theta)) nb_model$theta else 1
    start <- c(start, ln_alpha = log(1/theta))
    if (form == "nbp") {
      start <- c(start, P = 1.5)
    }
  } else {
    start <- start.vals
  }

  # --- Log-Likelihood Function ---
  log_likelihood <- function(params, y, X_fixed, X_rand, halton_draws, panel_index) {

    n_fixed <- ncol(X_fixed)
    n_rand <- ncol(X_rand)

    betas_fixed <- params[grep("^mean\\.", names(params))][1:n_fixed]
    betas_rand_mean <- params[grep("^mean\\.", names(params))][(n_fixed+1):(n_fixed+n_rand)]

    alpha <- exp(params["ln_alpha"])
    P <- if (form == "nbp") params["P"] else NULL

    mu_fixed <- as.vector(exp(X_fixed %*% betas_fixed))

    # Generate parameter draws for simulation
    draws <- matrix(NA, nrow = ndraws, ncol = n_rand)
    if (correlated) {
      chol_vals <- params[grep("^chol\\.", names(params))]
      chol_mat <- matrix(0, n_rand, n_rand)
      chol_mat[lower.tri(chol_mat, diag = TRUE)] <- chol_vals
      draws <- (qnorm(halton_draws) %*% t(chol_mat)) + matrix(betas_rand_mean, nrow = ndraws, ncol = n_rand, byrow = TRUE)
    } else {
      rand_sds <- abs(params[grep("^sd\\.", names(params))])
      for (i in 1:n_rand) {
        draws[, i] <- switch(rpardists[i],
                             "n"  = qnorm(halton_draws[, i], mean = betas_rand_mean[i], sd = rand_sds[i]),
                             "ln" = qlnorm(halton_draws[, i], meanlog = betas_rand_mean[i], sdlog = rand_sds[i]),
                             "t"  = qtri(halton_draws[, i], lower = betas_rand_mean[i] - rand_sds[i], upper = betas_rand_mean[i] + rand_sds[i], mode=betas_rand_mean[i]),
                             "u"  = qunif(halton_draws[, i], min = betas_rand_mean[i] - rand_sds[i], max = betas_rand_mean[i] + rand_sds[i]),
                             "g"  = qgamma(halton_draws[, i], shape = betas_rand_mean[i]^2 / rand_sds[i]^2, rate = betas_rand_mean[i] / rand_sds[i]^2)
        )
      }
    }

    # Calculate simulated probabilities
    xb_rand <- X_rand %*% t(draws)
    mu_sim <- mu_fixed * exp(xb_rand)

    size <- switch(form,
                   "nb1" = mu_sim / alpha,
                   "nb2" = 1 / alpha,
                   "nbp" = mu_sim^(2 - P) / alpha
    )

    # Probability for each observation and each draw
    prob_mat <- dnbinom(y, size = size, mu = mu_sim)
    prob_mat[prob_mat <= 0] <- 1e-300 # Avoid log(0)

    ## >> Improvement: Efficient panel likelihood calculation using base R
    # Sum of log-probabilities within each panel for each draw.
    # `rowsum` is highly optimized for this.
    panel_log_probs <- rowsum(log(prob_mat), group = panel_index, reorder = FALSE)

    # Average the panel likelihoods (not log-likelihoods) across draws
    # Use max of log-likelihoods for numerical stability (log-sum-exp trick)
    max_log_prob <- apply(panel_log_probs, 1, max)
    panel_probs <- exp(max_log_prob) * rowMeans(exp(panel_log_probs - max_log_prob))
    panel_probs[panel_probs <= 0] <- 1e-300

    log_panel_probs <- log(panel_probs)

    # For BHHH, maxLik needs one log-likelihood value per observation unit (panel)
    return(log_panel_probs)
  }

  # --- Model Fitting ---
  fit <- maxLik::maxLik(
    logLik = log_likelihood,
    start = start,
    method = method,
    control = list(iterlim = max.iters, printLevel = print.level),
    y = y,
    X_fixed = X_fixed,
    X_rand = X_rand,
    halton_draws = halton_draws,
    panel_index = panel_index
  )

  # --- Post-Estimation Processing ---
  result <- list()
  result$estimate <- fit$estimate
  result$vcov <- vcov(fit)
  result$logLik <- sum(logLik(fit)) ## >> Sum panel log-likelihoods for total
  result$nobs <- length(y)
  result$n_panels <- n_panels ## >> Store number of panels
  result$n_params <- length(start)

  result$formula <- formula
  result$rpar_formula <- rpar_formula
  result$panel_id <- panel ## >> Store panel id name
  result$form <- form
  result$correlated <- correlated
  result$rpardists <- rpardists
  result$fixed_names <- fixed_names
  result$rpar_names <- rpar_names

  if (correlated) {
    chol_vals <- fit$estimate[grep("^chol\\.", names(fit$estimate))]
    chol_mat <- matrix(0, length(rpar_names), length(rpar_names))
    chol_mat[lower.tri(chol_mat, diag = TRUE)] <- chol_vals
    result$rpar_vcov <- t(chol_mat) %*% chol_mat
    dimnames(result$rpar_vcov) <- list(rpar_names, rpar_names)
    result$rpar_cor <- cov2cor(result$rpar_vcov)
  }

  class(result) <- "rpnb"
  return(result)
}

# --- S3 Methods (Updated for Panel Info) ---

#' @export
print.rpnb <- function(x, ...) {
  cat("Random Parameters Negative Binomial Model\n")
  cat("-----------------------------------------\n")
  cat("NB Form:", toupper(x$form), "\n")
  cat("Correlated Parameters:", x$correlated, "\n")
  cat("Log-Likelihood:", round(x$logLik, 3), "\n")
  cat("Observations:", x$nobs, "\n")
  if (!is.null(x$panel_id)) {
    cat("Panels:", x$n_panels, "\n")
  }
  cat("\nCoefficients:\n")
  print(coef(x))
}

#' @export
summary.rpnb <- function(object, ...) {
  est <- object$estimate
  se <- sqrt(diag(object$vcov))
  z_val <- est / se
  p_val <- 2 * pnorm(-abs(z_val))

  results_df <- data.frame(
    Estimate = est,
    `Std. Error` = se,
    `z-value` = z_val,
    `Pr(>|z|)` = p_val,
    check.names = FALSE
  )

  cat("Random Parameters Negative Binomial Model\n")
  cat("-----------------------------------------\n")
  cat("NB Form:", toupper(object$form), "\n")
  cat("Log-Likelihood:", round(object$logLik, 3), "\n")

  ## >> Use number of panels for BIC if panel data exists
  n_obs_for_bic <- if (!is.null(object$panel_id)) object$n_panels else object$nobs
  cat("AIC:", round(-2 * object$logLik + 2 * object$n_params, 3), "\n")
  cat("BIC:", round(-2 * object$logLik + log(n_obs_for_bic) * object$n_params, 3), "\n\n")

  cat("Parameter Estimates:\n")
  print(results_df)

  if (object$correlated) {
    cat("\nRandom Parameter Correlation Matrix:\n")
    print(object$rpar_cor)
  }

  invisible(object)
}
