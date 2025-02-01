#' Logit
#'
#' @description
#' For a probability \eqn{p}, returns \eqn{log(p/(1-p))}.
#'
#'
#' @param p A probability \eqn{p}.
#'
#' @returns \eqn{log(p/(1-p))}
#' @export
#'
#' @examples
#'
#' p <- 0.34
#' logit(p)
#'
logit <- function(p){
  log(p / (1-p))
}
