#' Select variables based on posterior inclusion probabilities (PIPs)
#'
#' @param pip A vector or matrix of PIPs.
#' @param method Either \code{"fdr"} to perform Bayesian FDR correction (the
#' default) or \code{"cutoff"} to use a pre-defined threshold
#' @param fdr The FDR level when \code{method="fdr"}. Default is 0.05.
#' @param cutoff The PIP cutoff when \code{method="cutoff"}. Variables are
#' selected if their PIP is at least \code{cutoff}.
#'
#' @details
#' [metab_vs()] performs Bayesian variable selection via spike-and-slab priors of
#' the form \eqn{\beta_j=\tilde{\beta}_j\gamma_j}, \eqn{\tilde{\beta}\sim\mathcal{N}(0,\nu^2)},
#' \eqn{\gamma_j\in\{0,1\}}. Here, \eqn{\gamma_j} is an indicator variable
#' for whether \eqn{\beta_j} is zero. The posterior mean of \eqn{\gamma_j} is
#' called the posterior inclusion probability (PIP) and represents the conditional
#' probability that \eqn{X_{ji}} is associated with \eqn{V_i} given the data.
#' The purpose of this function is to classify variables as associated or
#' unassociated with the response variable by comparing their observed PIPs to
#' a cutoff value. A variable is declared to be "selected" if
#' its PIP is greater than or equal to the cutoff.
#'
#' The function can be applied in two ways: setting \code{method="cutoff"} to
#' choose the cutoff manually, or setting \code{method="fdr"} to calculate the
#' lowest cutoff that achieves the desired nominal false discovery rate (FDR).
#' The latter approach is based on the Bayesian FDR method proposed
#' by Newton et al. (2004, Biostatistics).
#'
#' @returns A logical object of the same dimensions as \code{pip} indicating
#' selected and non-selected variables.
#' @export
#'
#' @references
#' 1. Newton MA, Noueiry A, Deepayan S, Paul A. Detecting differential gene
#' expression with a semiparametric hierarchical mixture method. Biostatistics.
#' 2004;5(2):155-176.
#'

pip_to_variables <- function(pip, method = c("fdr", "cutoff"),
                             cutoff = 0.5, fdr = 0.05){

  method <- match.arg(method)
  pip_v <- c(pip) # if pip is a matrix, convert it to a vector

  if (all(pip_v == 0)) return(pip)

  # Determine cutoff
  if (method == "fdr"){
    cutoffs <- unique(pip_v) # potential cutoffs
    cutoffs <- cutoffs[cutoffs > 0]
    n_cutoffs <- length(cutoffs)
    n_variables <- length(pip_v)

    cutoff_mat <- matrix(cutoffs,
                         nrow = n_variables,
                         ncol = n_cutoffs,
                         byrow = TRUE)
    decisions <- pip_v >= cutoff_mat # selection indicator matrix for all pips and cutoffs
    fdrs <- colSums(decisions * (1 - pip_v)) / colSums(decisions) # cutoff-specific FDRs
    which_fdr <- which(fdrs <= fdr) # which cutoffs meet desired FDR?
    if (length(which_fdr) == 0){
      return(pip * 0)
    } else{
      cutoff <- min(cutoffs[which_fdr])
    }

  }

  (pip >= cutoff) * TRUE



}
