#' Standardize metabolites
#'
#' @param Y Numeric matrix of metabolite values. Each row is a subject and
#' each column is a metabolite. Values will be treated as point mass values
#' (PMVs) if they are \code{NA} or if they are below \code{psi}.
#' @param psi Numeric vector of metabolite-specific detection limits, length \code{ncol(Y)}.
#' Default is \code{apply(Y,2,min,na.rm=TRUE)}.
#' All \code{Y} values below their corresponding \code{psi} will be treated as
#' PMVs (in addition to any \code{Y} values that are \code{NA}).
#'
#' @return List of the following components:
#' * `Y`: A matrix of standardized metabolite values
#' * `means`: A vector of means used for standardization
#' * `sds`: A vector of standard deviations used for standardization.
#'
#' @details
#' Standardizing \code{Y} before analysis is essential for specifying hyperparameters.
#' However, this step is not trivial since \code{Y} is expected to contain
#' missing values called point mass values (PMVs) (see [metab_vs()]). To resolve
#' this, we impute \code{Y} using the half-minimum rule, standardize, then
#' re-introduce the missing values. This approach, though imperfect, is useful
#' for placing metabolites on a comparable numerical scale.
#'
#' @importFrom stats sd
#' @export
#'

standardize_metab <- function(Y, psi = apply(Y,2,min,na.rm=TRUE)){

  q <- ncol(Y)
  Y_na <-
    sapply(
      seq_len(q),
      function(x){
        (Y[,x] < psi[x] | is.na(Y[,x]))
      }
    )
  Y_imp <-
    sapply(
      seq_len(q),
      function(x){
        ifelse(Y_na[,x], psi[x]/2, Y[,x])
      }
    )
  means <- colMeans(Y_imp)
  sds <- apply(Y_imp,2,sd)
  Y <- as.matrix(scale(Y_imp))
  Y[Y_na] <- NA
  list(
    Y = Y,
    means = means,
    sds = sds,
    psi = (psi - means) / sds
  )
}
