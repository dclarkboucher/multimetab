## code to prepare `metabolites` dataset goes here
set.seed(1)
n <- 500
beta0 <- 5
delta <- 2
sigma <- 1
beta <- rep(c(1,0), each = 5)
beta_c <- c(2,1,0.5)
psi <- 4
rho <- 0.9

p <- length(beta)
s <- length(beta_c)
R <- abs(stats::cov2cor(LaplacesDemon::rinvwishart(p * 1.5, diag(p))))
diag(R) <- 0
X <- as.matrix(scale(matrix(rnorm(n * p), nrow = n, ncol = p)))
C <- as.matrix(scale(matrix(rnorm(n * s), nrow = n, ncol = s)))


y_initial <-
  beta0 + c(X %*% beta) + rnorm(n, sd = sigma) + c(C %*% beta_c) +
  abs(rnorm(n, sd = delta))

u <- rbinom(n, 1, rho)
y <- ifelse(y_initial < psi | u == 0, NA, y_initial)

metabolites <-
  list(
    y = y,
    X = X,
    C = C,
    R = R,
    psi = 4
  )

usethis::use_data(metabolites, overwrite = TRUE)
