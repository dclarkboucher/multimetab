#' Select 'eta' hyper-parameter for Markov random field prior
#'
#' @param y Numeric vector containing observations of the outcome variable.
#' @param X Numeric matrix containing observations of the predictor variables.
#' @param psi Numeric detection limit (optional). If specified, y values below this
#' value will be treated as left-censored values from the skewed normal distribution.
#' @param R Relationship matrix for Markov Random Field variable selection prior (optional).
#' If specified, must be a symmetric non-negative matrix with diagonal zero and
#' dimension \code{ncol(X)}. If not specified, the
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
#' @param eta_tol Error bound for eta. Default is \code{0.1}.
#' @param percentile_check Percentile used to check for phase transition. The default
#' is \code{0.1}, which means phase transition is determined
#' by tracking the 10th percentile of the \eqn{\hat{\gamma}} vector.
#' @param quantile_change_max Maximum allowable change in the quantile corresponding
#' to \code{percentile_check}. Default is 0.05. The algorithm identifies phase
#' transition if the \code{percentile_check} percentile of the \eqn{\hat{\gamma}}
#' vector is mode than \code{quantile_change_max} greater than the model with \eqn{eta=0}.
#'
#' @returns The largest \eqn{eta} value such that the model does not experience
#' phase transition.
#' @export
#'
#'
get_eta_pt <-
  function(y, X, psi,
           R = matrix(0,ncol(X), ncol(X)),
           hyper_params = NULL, burnin, draws, thinning = 5,
           refine_betas = TRUE, model_skewness = TRUE,
           adaptive = FALSE,
           adapt_prop = 0.25, theta2 = 2, vs = TRUE, pmax = ncol(X),
           pmax_draws = 10,
           eta_tol = 0.1,
           percentile_check = 0.10,
           quantile_change_max = 0.05
  ){

    n <- length(y)
    p <- ncol(X)
    if (quantile_change_max <= 0 | quantile_change_max >= 1){
      stop("'quantile_change_max' should be in (0,1).")
    }
    if (percentile_check <= 0 | percentile_check >= 1){
      stop("'percentile_check' should be in (0,1).")
    }
    if (eta_tol <=0) stop("'eta_tol' should be a small positive number")

    if (nrow(X) != n) stop ("Incompatible dimensions")
    if (p != nrow(R) | p != ncol(R)) stop("Incompatible dimensions")

    if (!is.null(hyper_params$eta)){
      if (hyper_params$eta != 0) stop("`hyper_params$eta` should be zero")
    }

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
          eta = 0,
          rho_0 = 2 * 0.9,
          rho_1 = 2 * 0.1

        )


    }

    # need to make hyper_parameters.

    # Start with eta = 0
    out0 <-
      metab_vs(y = y, X = X, psi = psi,
               R = R, hyper_params = hyper_params,
               burnin = burnin, draws = draws,
               thinning = 1, refine_betas = refine_betas,
               model_skewness = model_skewness,  adaptive = adaptive,
               adapt_prop = adapt_prop, theta2 = theta2,
               vs = vs, pmax = pmax,
               pmax_draws = pmax_draws)
    stat0 <- quantile(colMeans(out0$gamma), percentile_check)

    # Binary search for eta
    min <- 1e-4
    max <- 4
    nit <- ceiling(log((max - min) / eta_tol) / log(2))

    for (i in 1:nit){
      eta_temp <- mean(c(min,max))
      hyper_params$eta <- eta_temp
      out <-
        metab_vs(y = y, X = X, psi = psi,
                 R = R, hyper_params = hyper_params,
                 burnin = burnin, draws = draws,
                 thinning = thinning, refine_betas = refine_betas,
                 model_skewness = model_skewness,  adaptive = adaptive,
                 adapt_prop = adapt_prop, theta2 = theta2,
                 vs = vs, pmax = pmax,
                 pmax_draws = pmax_draws)

      if (!out$phase_transition){
        stat <- quantile(colMeans(out$gamma), percentile_check)
        if (stat - stat0 > 0.05){
          max <- eta_temp
        } else{
          min <- eta_temp
        }
      } else{
        max <- eta_temp
      }

    }

    return(min)

}
