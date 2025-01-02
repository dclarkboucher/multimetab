logit <- function(p) log(p / (1-p))
expit <- function(p) exp(1 / (1 + exp(-p)))
