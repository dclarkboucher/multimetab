#' Determine PIPs
#'
#' @param pip The pips
#' @param method the method
#' @param cutoff cutoff
#' @param fdr fdr
#'
#' @export
#'

pip_to_variables <- function(pip, method = c("fdr", "cutoff"),
                             cutoff = 0.5, fdr = 0.05){

  method <- match.arg(method)
  pip_v <- c(pip)

  if (all(pip_v == 0)) return(pip)

  # Determine cutoff
  if (method == "fdr"){
    cutoffs <- unique(pip_v)
    cutoffs <- cutoffs[cutoffs > 0]
    n_cutoffs <- length(cutoffs)
    n_variables <- length(pip_v)

    cutoff_mat <- matrix(cutoffs,
                         nrow = n_variables,
                         ncol = n_cutoffs,
                         byrow = TRUE)
    decisions <- pip_v >= cutoff_mat
    fdrs <- colSums(decisions * (1 - pip_v)) / colSums(decisions)
    which_fdr <- which(fdrs <= fdr)
    if (length(which_fdr) == 0){
      return(pip * 0)
    } else{
      cutoff <- min(cutoffs[which_fdr])
    }

  }

  (pip >= cutoff) * TRUE



}
