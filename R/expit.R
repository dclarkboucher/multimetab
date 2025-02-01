#' Expit
#'
#' Inverse of the \code{\link{logit}} function, defined as
#' \eqn{expit(x) = 1/(1 + exp(-x))}
#'
#' @param x A real number \eqn{x}
#'
#' @returns \eqn{1/(1 + exp(-x))}
#' @export
#'
#' @examples
#' p <- 0.3
#'
#' x <- logit(p)
#' expit(x)
#'
#'
expit <- function(x){

  1 / (1 + exp(-x))
}
