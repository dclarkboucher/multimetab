#' Bayesian variable selection for metabolite endpoints
#'
#' @param y Numeric vector of metabolite measurements. Values will be
#' treated as point mass values (PMVs) if they are \code{NA} or if they are below \code{psi}.
#' Can be "standardized" using [standardize_metab()] for ease of prior specification.
#' @param X Numeric matrix of the predictor variables that will
#' be considered for variable selection. Rows are the subjects and columns
#' are the variables. It is recommended to standard continuous variables to have
#' mean 0 and variance 1.
#' @param C Numeric matrix of confounding variables that are excluded
#' from variable selection (optional). Rows are the the subjects and columns are the variables.
#' It is recommended to standard continuous variables to have
#' mean 0 and variance 1. If provided, should NOT include an intercept term.
#' @param psi Numeric value representing the metabolite's detection limit. All
#' \code{y} values below \code{psi} will be treated as PMVs and inferred via data
#' augmentation. The default is \code{min(y, na.rm =TRUE)}, in which case only
#' \code{NA} values of \code{y} will be treated as PMVs.
#' @param R Relationship matrix for Markov Random Field variable selection prior (optional).
#' If specified, must be a symmetric non-negative matrix with diagonal zero and
#' dimension \code{ncol(X)}. If not specified, the
#' variable selection hyperprior will be independent and identically across
#' predictors, with prior inclusion probability \code{expit(omega)}. See details.
#' @param hyper_params List of hyper-parameters for the MCMC. See details.
#' @param burnin Number of burnin draws for the MCMC.
#' @param draws Number of post-burnin draws for the MCMC.
#' @param thinning Thinning parameter for the MCMC. Default is 5, meaning one in
#' every 5 posterior MCMC draws will be retained. Higher values reduce the correlation
#' between retained draws.
#' @param nchains Number of MCMC chains. Default is 1.
#' @param type Sampling scheme used for variable selection. Either \code{"gibbs"}
#' for Gibbs sampling (the default) or \code{"metro"} for Metropolis-Hastings
#' @param psamp If using Gibbs sampling, the number of selection indicators to update
#' in each iteration. 1 in 5 variables are sampled by default.
#' @param refine_betas Logical indicator for whether the regression coefficients
#' from the Metropolis-Hastings algorithm should be re-sampled with an additional
#' Gibbs step to improve mixing. Default is \code{TRUE}. If \code{FALSE} (not
#' recommended) and \code{type="metro"}, the regression parameters will be updated
#' using traditional stochastic search variable selection. If \code{type="metro"},
#' ignored.
#' @param theta2 Sampling variance of the proposed beta in the Metropolis-Hastings.
#' If \code{adaptive = TRUE}, this parameter is chosen adaptively for each parameter
#' based on the last \code{adapt_prop} proportion of the burnin draws. Ignored
#' unless \code{type="metro"}.
#' @param adaptive Logical value indicating whether to the Metropolis-Hastings
#' parameter \code{theta2} should be tuned based on samples drawn during the
#' burnin stage of the MCMC. Ignored unless \code{type="metro"}.
#' @param adapt_prop The proportion of burnin draws used to tune the Metropolis-Hastings
#' parameter. Ignored unless \code{adaptive=TRUE}.
#' @param vs Logical value indicating whether to perform variable selection. Default is \code{TRUE}.
#' @param zi Logical value indicating whether to fit a zero-inflated model. Default is \code{TRUE}.
#' @param model_skewness Logical value indicating whether to infer the
#' skewness parameter. Default is \code{TRUE}. If \code{FALSE}, the parameter
#' is fixed at zero and the error term is assumed to be normally distributed.
#' @param pmax Upper limit for the number of selected variables. If more than \code{pmax}
#' variables are selected in \code{pmax_draws} consecutive draws, the algorithm
#' is terminated early and indicates that the phase transition
#' boundary has been crossed. This is helpful for checking whether the prior
#' specification is inappropriate without needing to run the MCMC scheme to
#' convergence.The default value is \code{ncol(X)}, meaning there is no early stopping.
#' @param pmax_draws Used jointly with \code{pmax} to prevent phase transition. Default is \code{10}.
#'
#'
#'
#' @details
#'
#' ### Model Specification
#' \code{metab_vs} is the workhorse function for performing variable selection
#' with metabolites as endpoints. Let \eqn{Y_i} be the observed value of the
#' metabolite (possibly censored), let \eqn{\tilde{Y}_i} be the true value of the
#' metabolite (unknown if \eqn{Y_i} is censored), let \eqn{U_i} be an indicator
#' variable for whether \eqn{\tilde{Y}_i} is zero, and let \eqn{V_i} be the true
#' value of the metabolite when \eqn{\tilde{Y}_i} is not zero. We assume
#' \eqn{U_i} follows a Bernoulli distribution with probability \eqn{\rho}, and
#' we model \eqn{V_i} as
#' \deqn{V_i=\beta_0+\sum_{j=1}^p X_{ji}\beta_j + \sum_{k=1}^s C_{ki}\alpha_k + \epsilon_i}
#' where \eqn{\epsilon_i} follows a skew-normal (SN) distribution with variance
#' parameter \eqn{\sigma^2} and skewness parameter \eqn{\delta}, using the
#' parameterization by Sahu et al. (2009). We assume that missing values
#' arise when \eqn{\tilde{Y}_i<\psi}, where \eqn{\psi} is a positive detection limit
#' treated as known. Our primary interest lies in the regression coefficients
#' \eqn{\beta_1,\dots,\beta_p}, which represent the associations between the covariates
#' \eqn{X_{1i},\dots, X_{pi}} and the non-zero metabolites. We perform variable
#' selection among \eqn{X_{1i},\dots, X_{pi}} to determine which variables
#' are associated with the metabolite. The framework also allows
#' specifying a low-dimensional vector of confounding variables, \eqn{C_{1i},\dots,C_{si}},
#' that are fixed in the model and do not undergo variable selection.
#'
#' ### Prior distributions
#' The proposed framework performs variable selection via spike-and-slab priors of
#' the form \eqn{\beta_j=\tilde{\beta}_j\gamma_j}, \eqn{\tilde{\beta}\sim\mathcal{N}(0,\nu^2)},
#' \eqn{\gamma_j\in\{0,1\}}. Here, \eqn{\gamma_j} is an indicator variable
#' for whether \eqn{\beta_j} is zero. The posterior mean of \eqn{\gamma_j} is
#' called the posterior inclusion probability (PIP) and represents the conditional
#' probability that \eqn{X_{ji}} is associated with \eqn{V_i} given the data.
#'
#' The framework optionally places a sophisticated multivariate prior
#' distribution, the Markov random field (MRF) prior, on the vector
#' \eqn{\gamma=(\gamma_1,\dots,\gamma_p)^T}. Let \eqn{R} (\code{R}) denote
#' a non-negative, symmetric \eqn{p\times p} matrix with diagonal zero. The
#' prior has the form
#' \deqn{P(\gamma)\propto\exp(\omega\sum_j\gamma_j + \eta\gamma^TR\gamma)}
#' and is designed such that greater values of \eqn{R_{j,j^\prime}} result in greater
#' probability of selecting variables \eqn{j} and \eqn{j^\prime} jointly. This
#' is useful when there is prior belief that if one variable is associated with
#' the endpoint, then similar variables may be as well.
#' Within this prior, the hyper-parameter \eqn{\omega} controls the baseline
#' selection probability across all variables. The hyper-parameter \eqn{\omega}
#' controls the degree to which \eqn{R} informs variable selection. Whereas
#' \eqn{\omega} may be chosen based on prior belief about the proportion
#' of predictive variables (for example, setting \code{omega=logit(0.05)}),
#' \eqn{\eta} must be chosen carefully to avoid "phase transition behavior", a
#' phenomenon in which large \eqn{\eta} cause the model to select almost every
#' variable. We recommend choosing \eqn{\eta} using the function [get_eta_simple()],
#' which chooses a reasonable value of \eqn{\eta} based on \eqn{\omega} and
#' \eqn{R} and requires little computational burden. Alternatively, one may
#' use [get_eta_pt()] to approximate the lowest \eqn{\eta} value below the
#' phase transition point; however, this approach is computationally costly
#' since it requires sampling from the posterior distribution.
#' If \eqn{R} is not specified or if \eqn{\eta=0}, the MRF prior reduces to a
#' traditional independent hyper-prior in which the prior chance that \eqn{\gamma_j=1}
#' is \eqn{1/(1+\exp(-\omega))} (\code{expit(omega)}).
#'
#' For the remaining parameters, we use the standard conjugate priors:
#' \eqn{\beta_0\sim \mathcal{N}(0,\nu_0^2)}, \eqn{\sigma^2\sim \mathcal{IG}(0.5\xi_0,0.5\xi_0\sigma^2_0)},
#' \eqn{\delta\sim\mathcal{N}(0,\nu_d^2)}, \eqn{\alpha_k\sim\mathcal{N}(0,\lambda_k^2)}, and
#' \eqn{\rho\sim\mathcal{Beta}(\rho_0,\rho_1)}.
#'
#'
#'
#'
#' @return List of the following components:
#' * `samples`: A three-dimensional array of posterior samples organized by draw
#' index, chain, and parameter.
#' * `beta_hat`: A vector of the posterior means of the regression coefficients conditional
#' on their inclusion in the model, averaged over the chains. Variables that were
#' never selected will have a value of zero.
#' * `pip`: A vector of the posterior inclusion inclusion probabilities, averaged
#' over the chains.
#' * `estimates`: A matrix of posterior means of each parameter in each chain.
#' * `phase_transition`: An indicator of whether the phase transition boundary was
#' crossed. Will only equal \code{1} if \code{pmax} was specified and the algorithm
#' terminated early.
#'
#'
#' @useDynLib multimetab
#' @importFrom Rcpp evalCpp
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
  function(y, X, C = NULL, psi = min(y, na.rm = TRUE),
           R = matrix(0,ncol(X), ncol(X)),
           hyper_params = NULL, burnin, draws, thinning = 1, nchains = 1,
           type = c("gibbs", "metro"),
           psamp = ncol(X) %/% 5,
           refine_betas = TRUE,
           adaptive = FALSE,
           adapt_prop = 0.25, theta2 = 2,
           vs = TRUE, zi = TRUE, model_skewness = TRUE, pmax = ncol(X),
           pmax_draws = 10){

    # should maybe make theta2 a hyper-parameter
    # OR have one function used to create mcmc settings and one
    # function to create hyper-parameters. "control" like function, e.g.
    # Also, should this function determine eta by itself?
    # Maybe get_hyper_params function should have options for how to determine
    # eta. But then you would have to supply all the data.

    if (is.null(psi)) psi <- min(y, na.rm = TRUE)
    y[is.na(y)] <- psi - 1 # what matters is that the new values are below psi

    n <- length(y)
    p <- ncol(X)
    if (!is.null(C)){
      s <- ncol(C)
      XC <- cbind(X, C)
    } else{
      s <- 0
      XC <- X
    }

    type <- match.arg(type)
    gibbs <- type == "gibbs"

    if (is.null(hyper_params)){
      message("Using default hyper parameters")
      hyper_params <-
        list(
          nu2_0 = 10^2,
          nu2_d = 10^2,
          nu2 = rep(2^2, p),
          nu2_c = rep(10^2, s),
          xi_0 = 5,
          sigma2_0 = 4,
          omega = logit(0.05),
          eta = 0,
          rho_0 = 2 * 0.9,
          rho_1 = 2 * 0.1
        )
    }

    if ("lambda2" %in% names(hyper_params)){
      hyper_params$nu2_c = hyper_params$lambda2
    }

    if (nrow(X) != n) stop ("Incompatible dimensions")
    if (p != nrow(R) | p != ncol(R)) stop("Incompatible dimensions")

    total_draws = draws + burnin
    total_keep = floor(draws / thinning)

    # Output array
    nparams <- 2 * p + s + 4
    beta_names <- paste0("beta", 1:p)
    gamma_names <- paste0("gamma", 1:p)

    param_names <-
      c("beta0","sigma2","delta", "rho",
        beta_names, gamma_names)
    if (s > 0){
      beta_c_names <- paste0("beta_c", 1:s)
      param_names <- append(param_names, beta_c_names)
    }
    samples <- array(NA, dim = c(total_keep, nchains, nparams),
                dimnames = list(paste0("r", seq_len(total_keep)),
                                paste0("chain", seq_len(nchains)),
                                param_names
                ))

    # Table of posterior means
    estimates <- matrix(NA, nrow = nparams, ncol = nchains)
    colnames(estimates) <- paste0("chain", seq_len(nchains))
    rownames(estimates) <- param_names

    # MH moves
    moves_dat <- array(NA, dim = c(total_keep, nchains, 2),
                   dimnames = list(paste0("r", seq_len(total_keep)),
                                   paste0("chain", seq_len(nchains)),
                                   c("move","acceptance")))

    # Storage objects
    delta_samples <- numeric(total_keep)
    rho_samples <- numeric(total_keep)
    sigma2_samples <- numeric(total_keep)
    moves <- numeric(total_keep)
    acceptance <- numeric(total_keep)
    beta0_samples <- numeric(total_keep)
    gamma_samples <- matrix(0, total_keep, p)
    beta_samples <- matrix(0, total_keep, p)
    beta_c_samples <- matrix(0, total_keep, s)
    pt_check <- numeric(1)
    out <- list()

    for (j in seq_len(nchains)){
      bvs_mcmc(
        y = y, X = XC, s = s, psi = psi, hyper_params = hyper_params,
        R = R, theta2 = theta2, reps = draws, burnin = burnin,
        thinning = thinning, infer_delta = model_skewness,
        refine_betas = refine_betas, adaptive = adaptive, vs = vs, zi = zi,
        adapt_prop = adapt_prop, beta0_samples = beta0_samples,
        beta_samples = beta_samples, gamma_samples = gamma_samples,
        sigma2_samples = sigma2_samples, delta_samples = delta_samples,
        rho_samples = rho_samples, beta_c_samples = beta_c_samples,
        acceptance = acceptance, moves = moves,
        pt_check = pt_check, pmax = pmax, pmax_draws = pmax_draws,
        gibbs = gibbs, psamp = psamp
      )

      if (pt_check){
        message("Algorithm stopped early to prevent phase transition.")
        return(
          list(phase_transition = pt_check)
        )
      }

      samples[,j,"beta0"] <- beta0_samples
      samples[,j,"sigma2"] <- sigma2_samples
      samples[,j,"delta"] <- delta_samples
      samples[,j,"rho"] <- rho_samples
      samples[,j,gamma_names] <- gamma_samples
      samples[,j,beta_names] <- beta_samples
      if (s > 0){
        samples[,j,beta_c_names] <- beta_c_samples
      }
      estimates[,j] <- colMeans(samples[,j,])
      moves_dat[,j,"move"] <- moves
      moves_dat[,j,"acceptance"] <- acceptance

      # out[j] <-
      #   list(
      #     beta = beta_samples[,1:p],
      #     gamma = gamma_samples[,1:p],
      #     acceptance = acceptance,
      #     move = moves,
      #     beta0 = beta0_samples,
      #     delta = delta_samples,
      #     sigma2 = sigma2_samples,
      #     rho = rho_samples,
      #     phase_transition = pt_check
      #   )
      # if (s > 0){
      #   out$beta_c = beta_samples[,p + seq_len(s)]
      # }

    }

    dims <- ifelse(nchains > 1, 2, 1)
    pips <- colMeans(samples[,,gamma_names],
                     dims = dims)
    beta_mean <-
      colMeans( # Clever code to make beta NA if gamma = 0
        samples[,,gamma_names] / samples[,,gamma_names] * samples[,,beta_names],
        dims = dims, na.rm = TRUE
        )
    beta_mean[is.na(beta_mean)] <- 0
    names(beta_mean) <- beta_names
    # Generate means for each variable within chains
    output <-
    list(
      samples = samples,
      beta_hat = beta_mean,
      pip = pips,
      posterior_means = estimates,
      phase_transition = 0
    )
    if (!gibbs) output$moves <- moves_dat
    output
}
