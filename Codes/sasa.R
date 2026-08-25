
## Author: Ryan Martin (rgmartin@uic.edu)
## R codes to implement the SASA method for estimating finite mixtures.
## To run, first create the "sasa.so" file by typing, at the command line, 
##
##   R CMD SHLIB sasa.c
##
## With this file, the C functions "PR_finite", etc are accessible to R.
## See the corresponding paper for details of the method.


#===========================================================================================#

## Auxiliary functions

# Calculate PR estimate and likelihood for fixed support

prlik <- function(n, lik, f0, w, perm, nperm, penalty) {

  S <- length(f0)
  prlik <- penalty
  dyn.load("sasa.so")
  out <- .C("PR_finite", as.double(lik), as.integer(n), as.integer(S), as.double(w),
                         f=as.double(f0), lik=as.double(prlik), as.integer(nperm),
                         as.integer(perm))[c("f","lik")]
  dyn.unload("sasa.so")
  return(out)

}


# Calculate parametric likelihood function

likfn <- function(X, UU, d, ...) d(outer(X, 1 + 0 * UU), outer(1 + 0 * X, UU), ...)


# Proposal function for simulated annealing

sa.prop <- function(h, nbhd, r) {

  S <- length(h)
  if(sum(h) > 0) p <- 1 + h * (S / sum(h)) ** r else p <- 1 / S + numeric(S)
  k <- sample(nbhd, 1)
  iu <- sample(S, size=k, prob=p, replace=FALSE)
  h[iu] <- (h[iu] + 1) %% 2
  return(h)

}


# Function to count the modes of a kernel density estimate

mode.count <- function(x, adjust=1) {

    kd <- density(x, bw="nrd", adjust=adjust)
    diff <- as.numeric(kd$y[-1] - kd$y[-length(kd$y)] < 0)
    ddiff <- diff[-1] - diff[-length(diff)]
    count <- sum(ddiff == 1)

    return(count)

}



#===========================================================================================#


## MAIN FUNCTIONS

# Case 1: Kernel is fully known, e.g., Poisson or Gaussian with known scale.

sasa.known <- function(X, UU, start, d, gam, perm, nperm, nbhd, usereg, prob, x, sann.iter, ...) {

  n <- length(X)
  S <- length(UU)
  w <- 1 / (1:n + 1)**gam

  if(is.null(perm)) {

    perm <- c()
    for(j in 1:nperm) perm <- c(perm, sample(n) - 1)

  }

  lik <- likfn(X, UU, d, ...)

  L <- function(h) {

    if(all(h == 0)) o <- Inf else {

      pen <- usereg * (sum(h) * log(prob) + (S - sum(h)) * log(1-prob))
      lik0 <- lik[, h > 0]
      h0 <- h[h > 0]
      o <- prlik(n, lik0, h0 / sum(h0), w, perm, nperm, pen)$lik

    }

    return(o)

  }

  prop <- function(h) sa.prop(h, nbhd, 1)
  out <- optim(par=start, fn=L, gr=prop, method="SANN", control=list(maxit=sann.iter))
  h0 <- out$par
  pr.out <- prlik(n, lik[,h0 > 0], h0[h0 > 0] / sum(h0), w, perm, nperm, 0)
  fu <- pr.out$f
  if(is.null(x)) mx <- NULL
  else mx <- as.numeric(likfn(x, UU[h0 > 0], d, ...) %*% fu)

  return(list(U=UU[h0 > 0], fu=fu, x=x, mx=mx, nclusters=sum(h0), lik=pr.out$lik))

}


# General location-scale mixture

sasa.lsmix <- function(X, loc, scale, start, d, gam, nperm, perm, prob, x, sann.iter) {

  n <- length(X)
  N.loc <- length(loc)
  N.scale <- length(scale)
  w <- 1 / (1:n + 1) ** gam
  if(is.null(perm)) {

    perm <- c()
    for(j in 1:nperm) perm <- c(perm, sample(n) - 1)

  }
  DD <- outer(X, rep(1, N.loc)) - outer(rep(1, n), loc)

  L <- function(h) {

    if(sum(h > 0) == 0) oo <- Inf else {

      lik <- d(DD[,h > 0], 0, outer(rep(1, n), scale[h[h > 0]]))
      pen <- (sum(h > 0) * log(prob) + (N.loc - sum(h > 0)) * log(1-prob))
      oo <- prlik(n, c(lik), (h > 0)[h > 0] / sum(h > 0), w, perm, nperm, pen)$lik

    }

    return(oo)

  }

  prop <- function(h) {

    r <- 1
    if(all(h == 0)) p <- rep(1, N.loc)
    else p <- 1 + (h > 0) * ( 1 / mean(h > 0) )^r
    s <- sample(N.loc, size=1, prob=p, replace=FALSE)
    if(h[s] == 0) h[s] <- sample(N.scale, size=1)
    else if(runif(1) < sum(h > 0) / N.loc) h[s] <- 0
    else if(h[s] == 1) h[s] <- 2
    else if(h[s] == N.scale) h[s] <- N.scale - 1
    else h[s] <- sample(c(h[s]-1, h[s]+1), 1)
    return(h)

  }

  out <- optim(par=start, fn=L, gr=prop, method="SANN", control=list(maxit=sann.iter))
  h0 <- out$par
  loc0 <- loc[h0 > 0]
  scale0 <- scale[h0[h0 > 0]]

  lik <- d(outer(X, 1 + 0 * loc0), outer(rep(1, n), loc0), outer(rep(1, n), scale0))
  pr.out <- prlik(n, c(lik), (h0 > 0)[h0 > 0] / sum(h0 > 0), w, perm, nperm, 0)
  f0 <- pr.out$f

  if(is.null(x)) mx <- NULL
  else mx <- as.numeric(d(outer(x, 1 + 0 * f0), outer(1 + 0 * x, loc0), outer(1 + 0 * x, scale0)) %*% f0)

  return(list(loc=loc0, scale=scale0, f=f0, x=x, mx=mx, nclusters=sum(h0 > 0), lik=pr.out$lik))

}


#===========================================================================================#

## Galaxy data examples

# SASA estimation with Gaussian scale is fixed at sigma = 1.

galaxy.known <- function(sig=1, nperm=100, iter=5000, plot.it=FALSE) {

  X <- scan("galaxy.txt", quiet=TRUE)
  if(plot.it) x <- seq(min(X) - sd(X), max(X) + sd(X), len=300) else x <- NULL
  UU <- seq(5, 40, by=0.5)
  S <- length(UU)
  start <- 1 + 0 * UU
  np <- mode.count(X)
  out <- sasa.known(X, UU, start, dnorm, 1, NULL, nperm, 1, TRUE, np / S, x, iter, sd=sig)
  fu <- out$fu

  if(plot.it) {

    op <- par(mfrow=c(1,2), pty="m", mar=c(4.2, 4.2, 1, 1))
    plot(out$U, fu, type="h", lwd=2, xlab="u", ylab="f(u)", main="", xlim=c(5,40))
    hist(X, breaks=25, freq=FALSE, xlim=range(UU), col="gray", border="white", main="", ylab="m(y)", xlab="y")
    lines(out$x, out$mx, lwd=2)
    par(op)

  }

  return(out)

}


# SASA estimation of a Gaussian location-scale mixture

galaxy.glsmix <- function(nperm=100, iter=5000, plot.it=FALSE) {

  X <- scan("galaxy.txt", quiet=TRUE)
  if(plot.it) x <- seq(min(X) - sd(X), max(X) + sd(X), len=300) else x <- NULL
  loc <- seq(5, 40, by=0.5)
  scale <- seq(0.5, 1.5, by=0.1)
  start <- rep(1, length(loc))
  np <- max(1, mode.count(X))
  out <- sasa.lsmix(X, loc, scale, start, dnorm, 1, nperm, NULL, np / length(loc), x, iter)

  if(plot.it) {

    op <- par(mfrow=c(1,2), pty="m", mar=c(4.2, 4.2, 3, 1))
    plot(x=0, y=0, type="n", xlim=range(loc), ylim=range(scale), xlab=expression(mu), ylab=expression(sigma))
    abline(h=out$scale, v=out$loc, lty=3)
    text(x=out$loc, y=out$scale, paste(round(out$f, 2)), cex=0.8, font=2)
    hist(X, breaks=25, xlim=range(x), freq=FALSE, col="gray", border="white", main="", ylab="m(y)", xlab="y")
    lines(out$x, out$mx, lwd=2)

  }

  #print(cbind(loc=out$loc, scale=out$scale, f=out$f))
  return(out)

}
