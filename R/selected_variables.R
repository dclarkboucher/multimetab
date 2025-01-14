#' Identify selected variables from Bayesian model
#'
#'
#' @param fit A \code{list} as outputted by \code{metab_vs}.
#' @param method Either \code{"cutoff"} when providing a PIP cutoff to determine
#' significance, or \code{"fdr"} to determine the cutoff automatically using Bayesian FDR rules.
#' Default is \code{"fdr"}.
#' @param cutoff If \code{method = "cutoff"}, the PIP cutoff for declaring
#' significance. Default is \code{0.5}.
#' @param fdr If \code{method = "cutoff"}, the nominal FDR rate. Default is \code{0.05}.
#'
#' @returns A list containing selection indicators (\code{variables}) and the
#' PIP cutoff (\code{cutoff}).
#' @export
#'
selected_variables <- function(fit, method = c("fdr", "cutoff"),
                               cutoff = 0.5,fdr = 0.05){

  gamma_hat <- colMeans(fit$gamma)
  method <- match.arg(method)
  p <- length(gamma_hat)

  if (max(gamma_hat) == 0){
    message("No variables chosen in any draw.")
    return(
      list(
        variables = rep(0, p),
        cutoff = 1
      )
    )
  }


  if (method == "fdr"){

    # Search through potential cutoffs to find the best
    cutoffs <- unique(gamma_hat)
    cutoffs <- cutoffs[cutoffs > 0]
    c_mat <- matrix(cutoffs, nrow = p, ncol = length(cutoffs), byrow = TRUE)
    D <- gamma_hat >= c_mat
    fdr_observed <- colSums((1-gamma_hat) * D) / colSums(D)
    which_cutoffs_ok <- which(fdr_observed <= fdr)
    if (length(which_cutoffs_ok) == 0){
      message("No PIP cutoffs satisfying FDR.")
      return(
        list(
          variables = rep(0, p),
          cutoff = 1
        )
      )

    }
    cutoff <- min(cutoffs[which_cutoffs_ok])

  }

  list(
    variables = ifelse(gamma_hat >= cutoff, 1, 0),
    cutoff = cutoff
  )

}
