#' Bayesian variable selection for metabolite outcomes
#'
#' @param y Numeric vector containing observations of the outcome variable.
#' @param X Numeric matrix containing observations of the predictor variables.
#' @param psi Numeric detection limit (optional). If specified, y values below this
#' value will be treated as left-censored values from the skewed normal distribution.
#' @param R Relationship matrix for Markov Random Field variable selection prior (optional).
#' If specified, must be a symmetric numeric matrix with diagonal zero. If not specified, the
#' variable selection prior will be independent across predictors. See details.
#' @param hyper_params List of hyper-parameters for the MCMC. See details.
#' @param burnin Number of burnin draws for the MCMC.
#' @param draws Number of post-burnin draws for the MCMC.
#' @param thinning Thinning parameter for the MCMC. Default is 1, meaning every
#' post-burnin draw is retained.
#' @param refine_betas Logical value indicating whether the beta parameters be updated
#' separately from the gamma parameters using a Gibbs sampler. Default is TRUE.
#' @param model_skewness Logical value indicating whether to use the skew-normal
#' distribution rather than the normal distribution. Default is TRUE.
#' @param theta2 Sampling variance of the proposed beta in the Metropolis-Hastings. If adaptive,
#' this parameter is updated based on the last adapt_prop of the burnin draws.
#' @param adaptive Logical value indicating whether to tune the Metropolis-Hastings
#' parameter based on samples drawn during the burnin stage of the MCMC. Default
#' is false because my code for this isn't working.
#' @param adapt_prop The proportion of burnin draws used to tune the Metropolis-Hastings
#' parameter.
#' @param vs Logical value indicating whether to perform variable selection. Default is true.
#' @param pmax Upper limit for the number of selected variables. If more than \code{pmax} variable are selected
#' in \code{pmax_draws} consecutive draws, the algorithm is terminated early to prevent phase transition. The default
#' value is \code{ncol(X)}, meaning there is no early stopping.
#' @param pmax_draws Used with \code{pmax} to prevent phase transition. Default is \code{10}.
#'
#' @return List containing the following components:
#'
#'
#'
#' @export
#'
#' @examples
#' # Generate Data
#' n <- 50
#' beta <- c(1,1,0,0)
#' p <- length(beta)
#' beta0 <- 0
#' delta <- 3
#' sigma <- 1
#' rho <- 0.9
#' psi <- 0.5
#' X <- as.matrix(scale(matrix(rnorm(200), nrow = n, ncol = 4)))
#' u <- rbinom(n, size = 1, prob = rho)
#' y <- beta0 + c(X %*% beta) + sigma * rnorm(n) + delta * abs(rnorm(n))
#' out <- metab_vs(
#'     y = y, X = X, psi = psi,
#'     burnin = 100,
#'     draws = 100
#'
#'   )
#'
#'
metab_vs <-
  function(y, X, psi,
           R = matrix(0,ncol(X), ncol(X)),
           hyper_params = NULL, burnin, draws, thinning = 1,
           refine_betas = TRUE, model_skewness = TRUE,
           adaptive = FALSE,
           adapt_prop = 0.25, theta2 = 2, vs = TRUE, pmax = ncol(X),
           pmax_draws = 10){


  if (is.null(hyper_params)){
    message("Using default hyper parameters")
    hyper_params <-
      list(
        nu2_0 = 5^2,
        nu2_d = 5^2,
        nu2 = rep(2^2, p),
        xi_0 = 5,
        sigma2_0 = 4,
        omega = logit(0.05),
        eta = 0

      )


  }

  n <- length(y)
  p <- ncol(X)
  if (nrow(X) != n) stop ("Incompatible dimensions")
  if (p != nrow(R) | p != ncol(R)) stop("Incompatible dimensions")

  total_draws = draws + burnin
  total_keep = floor(draws / thinning)

  # Storage objects
  delta_samples <- numeric(total_keep)
  rho_samples <- numeric(total_keep)
  sigma2_samples <- numeric(total_keep)
  moves <- numeric(total_keep)
  acceptance <- numeric(total_keep)
  beta0_samples <- numeric(total_keep)
  gamma_samples <- matrix(0, total_keep, p)
  beta_samples <- matrix(0, total_keep, p)
  pt_check <- numeric(1)

  bvs_mcmc(
    y = y, X = X, psi = psi, hyper_params = hyper_params,
    R = R, theta2 = theta2, reps = draws, burnin = burnin,
    thinning = thinning, infer_delta = model_skewness,
    refine_betas = refine_betas, adaptive = adaptive, vs = vs,
    adapt_prop = adapt_prop, beta0_samples = beta0_samples,
    beta_samples = beta_samples, gamma_samples = gamma_samples,
    sigma2_samples = sigma2_samples, delta_samples = delta_samples,
    rho_samples = rho_samples, acceptance = acceptance, moves = moves,
    pt_check, pmax = pmax, pmax_draws = pmax_draws
  )

  return(
    list(
      beta = beta_samples,
      gamma = gamma_samples,
      acceptance = acceptance,
      move = moves,
      beta0 = beta0_samples,
      delta = delta_samples,
      sigma2 = sigma2_samples,
      rho = rho_samples,
      phase_transition = pt_check

    )
  )
}
