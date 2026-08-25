
##----------------------------------------------------------------------------------
## R codes for implementation of the predictive recursion (PR) method
## Author: Ryan Martin (rgmartin@uic.edu)
##
## Contents:
## * R+C version of PR
## * R-only version of PR
##
## The R+C versions, which run a bit faster than the R-only versions, require
## creation of a file "pr.so" containing the C routines to be called by R.  This
## is done by typing, at the shell prompt,
##
##   R CMD SHLIB pr.c
##
## In semiparametric mixture models, where additional non-mixing parameter are
## present, a PR marginal likelihood function can be used to estimate this non-
## mixing parameter.  For this, a PRML function can be defined and optimized
## using any standard routine, such as optim().  A general code for this can be
## written, but it's easy and probably better to do this case-by-case.
##
##----------------------------------------------------------------------------------


## Auxiliary numerical integration function

int <- function(f, x, tol=1e-10) {

  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) out <- sum((f[-1] + f[-n]) * (x[-1] - x[-n])) / 2
  else out <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n-1)]) + 2 * sum(f[seq(2, n-1, 2)]) )
  return(out)

}


## PR implemented with R and C
## (only implements continuous mixture; discrete mixture in "sasa.R" on website)

# INPUT:
# X = data
# d = parametric kernel, d(x, u)
# U = mixing density support (must be equispaced and of odd length)
# f0 = initial guess
# w = function to define weight sequence
# nperm = number of permuations
# perm = matrix of permutations of 1:length(X)
# ... = additional parameters for the kernel, e.g., a normal sd

# OUTPUT:
# f = PR estimate of f
# L = log-likelihood from PR fit

pr <- function(X, d, U, f0, w, nperm=1, perm=NULL, ...) {

  n <- length(X)
  ngrid <- length(U)
  if(ngrid %% 2 == 0) stop("grid g must have odd length")
  h <- U[2] - U[1]
  if(!all(U[-1] - U[-ngrid] == h)) stop("grid must be equispaced")
  if(missing(f0)) f0 <- 1 + 0 * U
  f0 <- f0 / int(f0, U)
  N <- nperm
  if(N == 1 && is.null(perm)) perm <- 1:n
  else if(N > 1 && is.null(perm)) {

    perm <- 1:n
    for(j in 1:(N-1)) perm <- c(perm, sample(n))

  }
  else perm <- c(perm)
  f.avg <- 0 * U
  L.avg <- 0
  if(missing(w)) w <- function(i) 1 / (i + 1)
  w <- w(1:n)
  lik <- c(d(outer(X, rep(1, ngrid)), outer(rep(1, n), U), ...))
  L <- rep(0, N)
  f <- rep(f0, N)
  dyn.load("pr.so")
  outpr <- .C("PR", as.double(lik), as.integer(n), as.integer(ngrid), as.double(h), as.double(w),
                    f=as.double(f), L=as.double(L), as.integer(N), as.integer(perm - 1))
  dyn.unload("pr.so")
  ftmp <- matrix(outpr$f, ngrid, N, byrow=FALSE)
  f.avg <- apply(ftmp, 1, mean)
  L.avg <- mean(outpr$L)
  return(list(f=f.avg, L=-L.avg))

}



## PR implemented with R only

# INPUT:
# X = data
# d = parametric kernel, d(x, u)
# U = mixing density support (must be equispaced and of odd length)
# f0 = initial guess
# w = function to define weight sequence
# nperm = number of permuations
# perm = matrix of permutations of 1:length(X)
# finite = logical, mixing dist is finite
# ... = additional parameters for the kernel, e.g., a normal sd

# OUTPUT:
# f = PR estimate of f
# L = log-likelihood from PR fit

pr.R <- function(X, d, U, f0, w, nperm=1, perm=NULL, finite=FALSE, ...) {

  n <- length(X)
  if(missing(f0)) f0 <- 1 + 0 * U
  if(missing(w)) w <- function(i) 1 / (i + 1)
  if(finite) f0 = f0 / sum(f0) else f0 = f0 / int(f0, U)
  N <- nperm
  if(N == 1 && is.null(perm)) perm <- 1:n
  else if(N > 1 && is.null(perm)) {

    perm <- matrix(0, n, N)
    perm[,1] <- 1:n
    for(j in 2:N) perm[,j] <- sample(n)

  }
  f.avg = 0 * U
  L.avg = 0
  for(j in 1:N) {

    f = f0
    L = 0
    x <- X[perm[,j]]
    for(i in 1:n) {

      num <- d(x[i], U, ...) * f
      if(finite) den <- sum(num) else den <- int(num, U)
      L <- L + log(den)
      f <- (1 - w(i)) * f + w(i) * num / den

    }

    f.avg <- (j - 1) * f.avg / j + f / j
    L.avg <- (j - 1) * L.avg / j + L / j

  }

  return(list(f=f.avg, L=-L.avg))

}

