
##-----------------------------------------------------------------------------------------
## TITLE:      imvch -- an inferential model (IM) approach for exact inference on
##             the heritiability coefficient in a linear mixed effect model.
##
## VERSION:    2nd version (07/05/2014).
##
## AUTHORS:    Q. Cheng, X. Gao, and R. Martin (rgmartin@uic.edu).
##
## REFERENCE:  Cheng, Gao & Martin (2014). "Exact prior-free probabilistic inference on
##             the heritability coefficient in a linear mixed model."  Electon. J. Stat.
##             (arXiv:1406.3521)
##-----------------------------------------------------------------------------------------

##~~~~~~~~~~ MAIN FUNCTION ~~~~~~~~~~##

# INPUT:
# eigen.out = eigenvalue output from 'vc.eigen' below
#     guess = value of rho to initialize MCMC (taken as the MLE in the example)
#     alpha = nominal confidence level
#    get.pi = compute the 100(1 - alpha)% plausibility interval?
#   plot.pl = compute the plausibility function of a grid and draw the plot?

# OUTPUT:
# list containing
# pl.int = the desired plausibility interval (NULL if get.pi == FALSE)
# pl.fun = plausibility values on a rho-grid (NULL if plot.pl == FALSE)

imvch <- function(eigen.out, guess, alpha=0.05, get.pi=TRUE, plot.pl=TRUE) {

  pl <- function(r) vch.pl(r, eigen.out, guess)
  max.R <- 1 - 1e-02
  RR <- seq(0, max.R, len=100)
  if(plot.pl) {

    o <- cbind(rho=RR, pl=sapply(RR, pl))
    plot(o[,1], o[,2], type="l", ylim=c(0, 1), xlab=expression(rho), ylab="Plausibility")

  } else o <- NULL
  if(get.pi) {

    f <- function(r) pl(r) - alpha
    if(f(0) > 0) lb <- 0 else {

      tmp <- guess
      while(f(tmp) < 0) tmp <- mean(c(tmp, 1))
      lb <- uniroot(f, c(0, tmp))$root

    }
    if(f(max.R) > 0) ub <- 1 else {

      tmp <- lb + 0.1
      while(f(tmp) < 0) tmp <- mean(c(tmp, lb))
      ub <- uniroot(f, c(tmp, max.R))$root

    }
    pl.int <- c(lb, ub)

  } else { pl.int <- NULL }
  return(list(pl.int=pl.int, out=o))

}


# Example 1 in the paper

assay.example <- function(alpha=0.05, get.pi=TRUE, plot.pl=TRUE) {

  load("assay.Rd")
  Y <- assay$Y
  X <- assay$X
  eigen.out <- vc.eigen(Y, X)
  mle <- 0.0306
  o <- imvch(eigen.out, mle, alpha, get.pi, plot.pl)
  print(o$pl.int)
  return(o)

}


##~~~~~~~~~~ AUXILIARY FUNCTIONS ~~~~~~~~~~##

# Computes the required eigenvalues, sufficient statistics, etc

vc.eigen <- function(y, X, Z, A) {

  if(!exists("ginv")) library(MASS)
  if(is.list(y)) {

    ng <- length(y)
    n <- sapply(y, length)
    N <- sum(n)
    Z <- matrix(0, nrow=N, ncol=ng)
    k <- 0
    for (g in 1:ng) {

      Z[k + (1:n[g]), g] <- 1
      k <- k + n[g]

    }
    Y <- unlist(y)
    if(missing(A)) A <- diag(ng)

  } else {

    Y <- y
    N <- length(Y)
    if(missing(Z)) stop("If y is not a list, then Z matrix is required")
    if(missing(A)) stop("If y is not a list, then A matrix is required")

  }
  if(missing(X)) X <- as.matrix(1 + numeric(N))
  p <- ncol(X)
  K <- svd(diag(N) - X %*% ginv(t(X) %*% X) %*% t(X), nu=N - p, nv=0)
  K <- t(K$u[, 1:(N - p)])
  KY <- K %*% Y
  G <- K %*% Z
  G <- G %*% A %*% t(G)
  e <- eigen(G)
  P <- e$vectors
  iValues <- round(e$values, 10)
  lambda <- unique(iValues)
  L <- length(lambda)
  S <- numeric(L)
  r <- integer(L)
  for (el in 1:L) {

    s <- which(iValues == lambda[el])
    r[el] <- length(s)
    a <- t(P[, s, drop=FALSE]) %*% KY
    S[el] <- sum(a * a)

  }
  return(list(lambda=lambda, L=L, r=r, S=S))

}


# Computes main ingredients for the IM for rho (Section 3.1 of the paper)

getDHT <- function(rho, eigen.out) {

  lambda <- eigen.out$lambda
  L <- eigen.out$L
  r <- eigen.out$r
  S <- eigen.out$S
  X <- (S[-L] / r[-L]) / (S[L] / r[L])
  f <- (1 + rho * (lambda[-L] - 1)) / (1 + rho * (lambda[L] - 1))
  g <- (lambda[-L] - lambda[L]) / (1 + rho * (lambda[-L] - 1)) / (1 + rho * (lambda[L] - 1))
  o.svd <- svd(matrix(g, ncol=L - 1), nv=L - 1)
  M <- t(o.svd$v[, -1, drop=FALSE])
  D <- rbind(M, 1 + 0 * g)
  T <- sum(log(X))
  H <- (D %*% log(X / f))[1:(L - 2)]
  return(list(D=D, Di=solve(D), H=H, T=T))

}


# Computes plausibility function via numerical integration of conditional density

vch.pl <- function(rho, eigen.out, guess) {

  lambda <- eigen.out$lambda
  L <- eigen.out$L
  r <- eigen.out$r
  DHT <- getDHT(rho, eigen.out)
  dtau <- function(tau) {

    z <- matrix(rep(DHT$H, length(tau)), nrow=L - 2, ncol=length(tau), byrow=FALSE)
    z <- drop(rbind(z, tau))
    z <- t(DHT$Di %*% z)
    rr <- r / 2
    lf <- z %*% (rr[-L] - 1) - sum(rr) * log(rr[L] + exp(z) %*% rr[-L]) + apply(z, 1, sum)
    return(exp(as.numeric(lf)))

  }
  f <- function(r) (1 + r * (lambda[-L] - 1)) / (1 + r * (lambda[L] - 1))
  v0 <- DHT$T - sum(log(f(guess)))
  p0 <- function(v) dtau(v) / dtau(v0)
  p <- function(v) p0(v) / integrate(p0, -Inf, Inf)$value
  vp <- function(v) v * p(v)
  mu <- integrate(vp, -Inf, Inf)$value
  P <- function(v) integrate(p, -Inf, v)$value
  phi <- sum(log(f(rho)))
  d.rho <- abs(DHT$T - phi - mu)
  pl <- 1 - (P(mu + d.rho) - P(mu - d.rho))
  return(pl)

}









