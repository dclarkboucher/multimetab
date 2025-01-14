logit <- function(p) log(p / (1-p))
expit <- function(x) 1 / (1 + exp(-x))
