#' Generate Various Types of Random and Quasi-Random Draws
#'
#' A unified interface for creating different types of draws for simulation,
#' including standard and scrambled Halton sequences based on the methods
#' described by Kolenikov (2012).
#'
#' @param n_draws The number of draws (or points) to generate.
#' @param n_dim The number of dimensions for each draw.
#' @param type A character string specifying the type of draws. Options are:
#'   \itemize{
#'     \item `"standard-halton"`: Standard, unscrambled Halton sequence.
#'     \item `"scrambled-halton-sqrt"`: Halton sequence scrambled with the square-root method.
#'     \item `"scrambled-halton-neg-sqrt"`: Halton sequence scrambled with the negative square-root method.
#'     \item `"scrambled-halton-rand-mult"`: Halton sequence scrambled with a random multiplier.
#'     \item `"scrambled-halton-rand-perm"`: Halton sequence scrambled with a random permutation.
#'     \item `"scrambled-halton-atanassov"`: Halton sequence scrambled with Atanassov's method.
#'     \item `"standard-sobol"`: Standard Sobol sequence.
#'     \item `"scrambled-sobol"`: Scrambled Sobol sequence.
#'     \item `"pseudo-random"`: Standard pseudo-random draws from a uniform distribution.
#'   }
#' @param ... Additional arguments passed to specific draw functions (e.g., `seed`
#'   for scramblers that use randomization).
#'
#' @return A matrix of draws with `n_draws` rows and `n_dim` columns.
#'
#' @importFrom randtoolbox halton sobol
#' @importFrom stats runif
#'
#'
#' @examples
#' # Generate 100 draws in 6 dimensions (5 & 6 are correlated)
#'
#' # Standard Halton sequence
#' draws1 <- make_draws(100, 6, type = "standard-halton")
#'
#' # Scrambled Halton sequence using the square-root method
#' draws2 <- make_draws(100, 6, type = "scrambled-halton-neg-sqrt")
#' #'
#' # Scrambled Halton sequence with a random permutation and a seed
#' draws3 <- make_draws(100, 6, type = "scrambled-halton-rand-perm", seed = 123)
#'
#' # Scrambled Halton sequence with atanassov
#' draws4 <- make_draws(100, 6, type = "scrambled-halton-atanassov")
#'
#' # Plot to compare
#' par(mfrow = c(1, 4))
#' plot(draws1[,5:6], main = "Standard Halton")
#' plot(draws2[,5:6], main = "Scrambled (Neg-Sqrt)")
#' plot(draws3[,5:6], main = "Scrambled (Rand Perm)")
#' plot(draws4[,5:6], main = "Scrambled (Atanassov)")
#'
#' @export
make_draws <- function(n_draws, n_dim, type, ...) {
  draws <- switch(
    type,
    `standard-halton` = randtoolbox::halton(n = n_draws, dim = n_dim),
    `scrambled-halton-sqrt` = .scrambled_halton_engine(n_draws, n_dim, .scrambler_sqrt, ...),
    `scrambled-halton-neg-sqrt` = .scrambled_halton_engine(n_draws, n_dim, .scrambler_neg_sqrt, ...),
    `scrambled-halton-rand-mult` = .scrambled_halton_engine(n_draws, n_dim, .scrambler_random_multiplier, ...),
    `scrambled-halton-rand-perm` = .scrambled_halton_engine(n_draws, n_dim, .scrambler_random_permutation, ...),
    `scrambled-halton-atanassov` = .scrambled_halton_engine(n_draws, n_dim, .scrambler_atanassov, ...),
    `standard-sobol` = randtoolbox::sobol(n = n_draws, dim = n_dim, scrambling = 0, ...),
    `scrambled-sobol` = randtoolbox::sobol(n = n_draws, dim = n_dim, scrambling = 3, ...),
    `pseudo-random` = matrix(stats::runif(n_draws * n_dim), nrow = n_draws, ncol = n_dim),
    stop("Invalid 'type' specified. See ?make_draws for options.")
  )
  return(draws)
}


#-------------------------------------------------------------------------------
#
#                 INTERNAL HALTON SCRAMBLING ENGINE
#
#-------------------------------------------------------------------------------

#' Core Engine for Generating Scrambled Halton Sequences
#'
#' This internal function generates Halton sequences scrambled by a user-specified
#' function, based on the methods described by Kolenikov (2012).
#'
#' @param n_draws The number of points in the sequence.
#' @param n_dim The number of dimensions.
#' @param scrambler_fun A function that performs the scrambling.
#' @param ... Additional arguments passed on to the scrambler function.
#' @return A matrix of scrambled Halton draws.
#' @noRd
.scrambled_halton_engine <- function(n_draws, n_dim, scrambler_fun, ...) {
  primes <- .generate_primes(n_dim)
  h <- matrix(0, nrow = n_draws, ncol = n_dim)
  i <- 1:n_draws

  # Loop over each dimension (prime)
  for (d in 1:n_dim) {
    p <- primes[d]
    # Vectorized calculation of base-p digits for all numbers
    max_digits <- if (n_draws > 0) floor(log(max(i), base = p)) + 1 else 1
    p_powers <- p^(0:(max_digits - 1))

    # Correctly create the digits matrix using an outer product.
    # This avoids the vector recycling that caused the error.
    digits <- floor(outer(i - 1, 1 / p_powers)) %% p

    # Apply the user-specified scrambling function
    scrambled_digits <- scrambler_fun(digits = digits, base = p, ...)

    # Vectorized radical inverse calculation
    inv_p_powers <- 1 / (p * p_powers)
    h[, d] <- scrambled_digits %*% inv_p_powers
  }
  return(h)
}

#-------------------------------------------------------------------------------
#
#                 KOLENIKOV (2012) SCRAMBLING FUNCTIONS
#
#-------------------------------------------------------------------------------

#' @noRd
.scrambler_sqrt <- function(digits, base, ...) {
  multiplier <- floor(sqrt(base))
  return((digits * multiplier) %% base)
}

#' @noRd
.scrambler_neg_sqrt <- function(digits, base, ...) {
  multiplier <- base - round(sqrt(base))
  return((digits * multiplier) %% base)
}

#' @noRd
.scrambler_random_multiplier <- function(digits, base, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed + base) # Vary seed by base
  multiplier <- sample(1:(base - 1), 1)
  return((digits * multiplier) %% base)
}

#' @noRd
.scrambler_random_permutation <- function(digits, base, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed + base) # Vary seed by base
  perm <- c(0, sample(1:(base - 1)))
  return(matrix(perm[digits + 1], nrow = nrow(digits), ncol = ncol(digits)))
}

#' @noRd
.scrambler_atanassov <- function(digits, base, ...) {
  K <- matrix(c(
    2, 1, 3, 1, 5, 4, 7, 2, 11, 9, 13, 9, 17, 2, 19, 1, 23, 13,
    29, 6, 31, 22, 37, 7, 41, 37, 43, 36, 47, 36, 53, 39, 59, 4,
    61, 26, 67, 13, 71, 12, 73, 35, 79, 66, 83, 60, 89, 68,
    97, 63, 101, 47, 103, 15, 107, 104, 109, 4, 113, 64
  ), ncol = 2, byrow = TRUE)

  k_val <- K[K[, 1] == base, 2]
  if (length(k_val) == 0) k_val <- 1 # Default if not found

  max_digits <- ncol(digits)
  j <- 0:(max_digits - 1)
  multipliers <- (k_val^j) %% base

  scrambled <- (digits * matrix(multipliers, nrow = nrow(digits), ncol = max_digits, byrow = TRUE)) %% base
  return(scrambled)
}


#-------------------------------------------------------------------------------
#
#                 INTERNAL HELPER FUNCTIONS
#
#-------------------------------------------------------------------------------

#' Generate the first n prime numbers
#' @param n The number of primes to generate.
#' @return A numeric vector of the first n primes.
#' @noRd
.generate_primes <- function(n) {
  if (n < 1) return(numeric(0))
  primes <- numeric(n)
  primes[1] <- 2
  if (n == 1) return(primes)
  count <- 1
  num <- 3
  while (count < n) {
    is_prime <- TRUE
    limit <- sqrt(num)
    for (p in primes[1:count]) {
      if (p > limit) break
      if (num %% p == 0) {
        is_prime <- FALSE
        break
      }
    }
    if (is_prime) {
      count <- count + 1
      primes[count] <- num
    }
    num <- num + 2 # Check only odd numbers
  }
  return(primes)
}

