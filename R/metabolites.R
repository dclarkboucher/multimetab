#' Simulated metabolite data
#'
#' Metabolite data simulated from a skew-normal mixture model.
#'
#' @format ## `metabolites`
#' A list with two elements:
#' \describe{
#'   \item{y}{Metabolite values for 100 subjects}
#'   \item{X}{Numeric matrix of 10 predictor variables}
#'   \item{C}{Numeric matrix of 3 confounding variables}
#'   \item{psi}{Detection limit of the metabolite. \code{Y} values below
#'   \code{psi} are censored}
#'   \item{R}{Covariate relationship matrix used in Markov Random Field prior}
#'   ...
#' }
"metabolites"
