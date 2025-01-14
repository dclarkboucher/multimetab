## code to prepare `metabolites` dataset goes here
set.seed(1)
n <- 100
beta0 <- 5
delta <- 2
sigma <- 1
beta <- rep(c(1,0), each = 5)
psi <- 4
rho <- 0.9

p <- length(beta)
R <- stats::cov2cor(LaplacesDemon::rinvwishart(p * 1.5, diag(p)))
diag(R) <- 0
X <- as.matrix(scale(matrix(rnorm(n * p), nrow = n, ncol = p)))

y_initial <-
  beta0 + c(X %*% beta) + rnorm(n, sd = sigma) +
  abs(rnorm(n, sd = delta))

u <- rbinom(n, 1, rho)
y <- u * (y_initial > psi) * y_initial

metabolites <-
  list(
    y = y,
    X = X,
    R = R,
    psi = 4
  )

usethis::use_data(metabolites, overwrite = TRUE)
