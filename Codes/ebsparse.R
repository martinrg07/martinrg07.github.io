
##-----------------------------------------------------------------------------------------------
## TITLE:        ebsparse -- Empirical Bayes estimation of sparse normal mean.
##
## VERSION:      3rd version (09/24/2014).
##
## AUTHORS:      R. Martin (rgmartin@uic.edu) and S.G. Walker.
##
## REFERENCE:    Martin and Walker (2014). "Asymptotically minimax empirical Bayes estimation
##               of a sparse normal mean vector," Electron. J. Stat. (arXiv:1304.7366)
##
## DESCRIPTION:  Method runs a Gibbs sampler to simulate from the empirical Bayes posterior
##               distribution described in the paper.  From this posterior, one can readily
##               compute the empirical Bayes estimator of the mean vector or any other feature
##               of the posterior; for example, in the paper, the empirical Bayes posterior
##               distribution for an associated weight parameter is shown, as well as empirical
##               Bayes posterior posterior inclusion probabilities.
##
## NOTE:         Authors thank Roger Koenker for suggestions to improve Gibbs sampler code.
##-----------------------------------------------------------------------------------------------

# INPUTS:
# X = data vector, length n
# k = kappa, value between 0 and 1 (close to 1)
# a = alpha, value should be relatively small; if NULL, use a method of moments estimate
# v = sigma^2, should be > 1 / (1 - k) for proper prior, but smaller values work
# M = number of Monte Carlo samples

eb.gibbs <- function(X, k, a, v, M) {

  n <- length(X)
  theta <- matrix(0, M, n)
  W <- numeric(M)
  if(is.null(a)) {

    D.hat <- sum(abs(X) <= sqrt(2 * log(n)))
    a <- D.hat / n / (n - D.hat)

  }
  for(m in 1:M) {

    if(m == 1) D <- sum(abs(X) <= sqrt(2 * log(n))) else D <- sum(theta[m-1,] == 0)
    w <- rbeta(1, a * n + D, 1 + n - D)
    W[m] <- w
    u <- runif(n)
    p1 <- w * exp(-k * X**2 / 2)
    p2 <- (1 - w) / sqrt(1 + k * v)
    p <- p1 / (p1 + p2)
    theta[m,] <- (u > p) * rnorm(n, X, sqrt(v / (1 + k * v)))

  }
  return(list(W=W, theta=theta))

}


# Simulation study

mse.sim <- function(n, s, A, reps, k, a, v, M) {

  set.seed(7)
  mse <- 0 * outer(s, A)
  for(i in seq_along(s)) {

    for(j in seq_along(A)) {

      true <- c(rep(A[j], s[i]), rep(0, n-s[i]))
      for(r in 1:reps) {

        X <- true + rnorm(n)
        est <- apply(eb.gibbs(X, k, a, v, M)$theta, 2, mean)
        mse[i,j] <- mse[i,j] + sum((est - true)**2) / reps

      }

    }

  }
  rownames(mse) <- s
  colnames(mse) <- A
  return(round(mse))

}

# To replicate simulations in Table 1 of the paper:
# mse.sim(n=200, s=c(10, 20, 40), A=7:8, reps=100, k=0.99, a=0.25, v=100, M=1000)








