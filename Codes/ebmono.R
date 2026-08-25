
##-----------------------------------------------------------------------------------------------
## TITLE:        ebmono -- Empirical Bayes inference on a monotone non-increasing density
##
## VERSION:      1st version (07/30/2018).
##
## AUTHOR:       R. Martin (rgmarti3@ncsu.edu) 
##
## REFERENCE:    Martin, "Empirical priors and posterior concentration rates for a monotone 
##               density," arXiv:1706.08567
##
## DESCRIPTION:  This code implements two methods for inference on a monotone non-decreasing 
##               density.  One is the empirical Bayes method proposed in the above reference; 
##               the other is a Bayesian Dirichlet process mixture model.  Implementation of 
##               the former is based on my own Gibbs sampler formulation, while the latter is
##               based on the slice sampler in Kalli, Griffin, and Walker (2011).
##
## DEPENDS ON:   'grenander' function in package 'fdrtool'.
##
## CONTENTS:     Two main functions (ebmono & dpmono), an example, and some auxiliary functions.
##-----------------------------------------------------------------------------------------------


library("fdrtool")

## MAIN FUNCTION 1: empirical Bayes 
## Inputs:
## X = data
## M = number of MCMC runs (burn-in of 0.2 * M is added automatically)
## xgrid = grid of x values at which density is estimated
## cc and dd = hyperparameters described in the paper

ebmono <- function(X, M, xgrid, cc, dd) {
  
  n <- length(X)
  mle <- mixunif.mle(X)
  S <- mle$S
  mu.mle <- mle$mu
  omega.mle <- mle$omega
  f.mle <- dmixunif(xgrid, omega.mle, mu.mle)
  burn <- round(0.2 * M)
  MM <- M + burn
  ff <- matrix(0, nrow=M, ncol=length(xgrid))
  omega <- omega.mle
  mu <- mu.mle
  alpha <- 1 + cc * omega.mle
  Z <- numeric(n)
  N <- numeric(S)
  for(m in 1:MM) {
    
    p <- outer(1:n, 1:S, function(i, s) omega[s] * dunif(X[i], 0, mu[s]))
    for(i in 1:n) Z[i] <- sample(S, size=1, prob=p[i,])
    for(s in 1:S) {
      
      is <- which(Z == s)
      N[s] <- length(is)
      if(N[s] == 0) mu[s] <- rpareto(dd, mu.mle[s]) 
      else {
        
        max.s <- max(c(mu.mle[s], X[is]))
        mu[s] <- rpareto(N[s] + dd, max.s) 
        
      }
      
    }
    omega.tmp <- rgamma(S, N + alpha, 1)
    omega <- omega.tmp / sum(omega.tmp)
    if(m > burn) ff[m - burn,] <- dmixunif(xgrid, omega, mu)
    
  }
  f.post <- colMeans(ff)
  f.upper <- apply(ff, 2, quantile, probs=0.975)
  f.lower <- apply(ff, 2, quantile, probs=0.025)
  return(list(ff=ff, f.mle=f.mle, f.post=f.post, f.lower=f.lower, f.upper=f.upper))
  
}


## MAIN FUNCTION 2: Dirichlet process mixture model
## Inputs:
## X = data
## M = number of MCMC runs (burn-in of 0.2 * M is added automatically)
## xgrid = grid of x values at which density is estimated
## alpha = DP precision parameter
## shape = shape parameter on inverse gamma base measure
## scale = scale parameter on inverse gamma base measure
## kappa = input for the slice sampling algorithm 

dpmono <- function(X, M, xgrid, alpha=1, shape=NULL, scale=NULL, kappa=0.5) {
  
  if(any(is.null(c(shape, scale)))) {
    
    m.X <- mean(X)
    v.X <- var(X)
    a <- 2 + m.X**2 / v.X
    b <- (a - 1) * m.X
    
  } else { a <- shape; b <- scale }
  n <- length(X)
  B <- round(0.2 * M)
  fx <- 0 * xgrid
  d <- mu.tmp <- numeric(n)
  mu <- as.numeric(quantile(X, 1:5 / 5))
  for(i in 1:n) d[i] <- min(which(X[i] <= mu))
  g <- function(j) (1 - kappa) * kappa**(j - 1)
  ff <- matrix(0, nrow=M, ncol=length(xgrid))
  for(r in 1:(M + B)) {
    
    u <- runif(n, 0, g(d))
    L <- 1 + floor((log(u) - log(1 - kappa)) / log(kappa))
    J <- max(d)
    NN <- max(c(L, J))
    xi <- g(1:NN)
    v <- w <- numeric(NN)
    if(NN > J) mu <- c(mu, numeric(NN - J))
    for(j in 1:NN) {
      
      XX <- X[d == j]
      nj <- length(XX)
      if(nj > 0) mu[j] <- 1 / rtrgamma(a + nj, 1 / b, 1 / max(XX))
      else mu[j] <- 1 / rgamma(1, a, 1 / b)
      v[j] <- rbeta(1, 1 + sum(d == j), alpha + sum(d > j))
      if(j == 1) w[j] <- v[j] else w[j] <- v[j] * prod(1 - v[1:(j - 1)])
      if(r > B) ff[r - B,] <- ff[r - B,] + w[j] * dunif(xgrid, 0, mu[j])
      
    }
    for(i in 1:n) {
      
      p <- (xi > u[i]) * w * dunif(X[i], 0, mu) / xi
      mu.tmp[i] <- mu[sample(NN, size=1, prob=p)]
      
    }
    mu <- unique(mu.tmp)
    d <- apply(outer(mu.tmp, mu, "-"), 1, function(x) which(x == 0))
    
  }
  f.post <- colMeans(ff)
  f.upper <- apply(ff, 2, quantile, probs=0.975)
  f.lower <- apply(ff, 2, quantile, probs=0.025)
  return(list(ff=ff, f.post=f.post, f.lower=f.lower, f.upper=f.upper))
  
}



## EXAMPLE

if(FALSE) {
  
  n <- 100
  M <- 2000
  X <- rexp(n, 1)
  f.true <- function(x) dexp(x)
  xgrid <- seq(0, max(X), len=100)
  cc <- 0.01 * n**(5 / 3) / (log(n))**(2 / 3)
  dd <- log(n) / 20
  eb.time <- system.time(o <- ebmono(X, M, xgrid, cc, dd)); print(eb.time)
  ylim <- c(0, 1.3)
  plot(xgrid, 0 * xgrid, type="n", xlab="x", ylab="Density", ylim=ylim)
  for(m in 1:M) lines(xgrid, o$ff[m,], col="lightgray")
  lines(xgrid, o$f.mle, col="blue")
  lines(xgrid, f.true(xgrid), col="red")
  lines(xgrid, o$f.post)
  lines(xgrid, o$f.lower, lty=2)
  lines(xgrid, o$f.upper, lty=2)
  dp.time <- system.time(oo <- dpmono(X, M, xgrid)); print(dp.time)
  plot(xgrid, 0 * xgrid, type="n", xlab="x", ylab="Density", ylim=ylim)
  for(m in 1:M) lines(xgrid, oo$ff[m,], col="lightgray")
  lines(xgrid, o$f.mle, col="blue")
  lines(xgrid, f.true(xgrid), col="red")
  lines(xgrid, oo$f.post)
  lines(xgrid, oo$f.lower, lty=2)
  lines(xgrid, oo$f.upper, lty=2)
  
}


## AUXILIARY FUNCTIONS:
## The first extracts the mixture components from Grenander's estimate;
## the second evaluates a mixture of uniforms density for given mixture components;
## the third simulates from a Pareto distribution;
## the fourth simuluates from a truncated gamma using rejection sampling.

mixunif.mle <- function(X) {
  
  F <- ecdf(X)
  o <- grenander(F)
  mu <- o$x.knots
  f <- o$f.knots
  S <- length(mu)
  d <- sum(f[-S] * diff(mu))
  ff <- 0 * f
  ff[1] <- max((1 - d) / mu[1], f[1])
  ff[-1] <- f[-S]
  omega <- numeric(S)
  omega[S] <- ff[S] * mu[S]
  for(s in S:2) omega[s-1] = (ff[s-1] - ff[s]) * mu[s-1]
  return(list(S=S, omega=omega, mu=mu))
  
}


dmixunif <- function(x, omega, mu) {
  
  f <- function(z) sum(omega * dunif(z, 0, mu))
  return(sapply(x, f))
  
}


rpareto <- function(pow, cut) {
  
  u <- runif(1)
  return(cut / (1 - u)**(1 / pow))
  
}


rtrgamma <- function(shape, scale, upper) {
  
  a <- shape
  b <- scale
  z <- upper 
  Fz <- pgamma(z, shape=a, scale=b)
  if(Fz == 0) Y <- runif(1, 0.95 * z, z)
  else if(Fz > 1e-06) Y <- qgamma(Fz * runif(1), shape=a, scale=b) 
  else {
    
    r <- function(y) z * dgamma(y, shape=a, scale=b) * (y <= z) / Fz
    if(a < 1) stop("need gamma shape >= 1") else M <- r(min(b * (a - 1), z))
    repeat {
      
      YY <- runif(1, 0, z)
      if(runif(1) <= r(YY) / M) { Y <- YY; break }
      
    }
    
  }
  return(Y)
  
}



