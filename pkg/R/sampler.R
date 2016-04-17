bootstrap <- function(B = 500) {
  sampfun <- function(n) replicate(B, sample(1L:n, n, replace = TRUE))
  list(method = "Bootstrap sampling", sampler = sampfun)
}

subsampling <- function(B = 500, frac = 0.8) {
  sampfun <- function(n) replicate(B, sample(1L:n, floor(frac * n), replace = FALSE))
  list(method = paste0("Subsampling with ", sprintf("%1.1f", 100 * frac, 2), "% data"), 
    sampler = sampfun)
}

crossvalidation <- function(k = 5) {
  sampfun <- function(n) {
    ret <- matrix(NA, ncol = ceiling(n/k), nrow = k)
    ret[1L:n] <- sample(1L:n)
    t(ret)
  }
  list(method = paste0(k, "-fold cross validation"), sampler = sampfun)
}

jackknife <- function(d = 1, maxrep = 5000) {
  sampfun <- function(n) {
    if (choose(n, d) > maxrep) 
      stop("Maximum number of repetitions allowed exceeded. Reduce d!")
    apply(utils::combn(1L:n, d), 2, function(x) (1L:n)[-x])
  }
  list(method = paste0("Leave-", d, "-out jackknife"), sampler = sampfun)
} 
