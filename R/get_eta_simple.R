get_eta_simple <-
  function(omega, R, etas_try, burnin = 10000,
           draws = 30000, thinning = 5, quantile = 0.95){

    out <-
      get_prior_stat(omega, etas_try, R, draws = draws,
                     burnin = burnin, thin = thinning, quantile = quantile) |>
      as.data.frame()

    # Use out to determine prior.
    p <- ncol(R)
    target <- qbinom(quantile, p, prob = logit(2 * expit(omega)))
    eta_indices <- which(out$quantile < target)
    if (length(eta_indices) == 0) stop("No etas found. Try new eta_list.")
    max(etas_try[eta_indices])
}
