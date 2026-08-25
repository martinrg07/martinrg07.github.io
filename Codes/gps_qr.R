
##-----------------------------------------------------------------------------------------------
## TITLE:        gps.qr -- Scaling the Gibbs posterior in quantile regression via GPS.
##
## VERSION:      1st version (09/03/2015).
##
## AUTHORS:      N. Syring and R. Martin (rgmartin@uic.edu).
##
## REFERENCE:    Scaling the Gibbs posterior credible regions, arXiv:
##
## DESCRIPTION:  Employs a combination of bootstrap and stochastic approximation to scale the
##               Gibbs posterior in such a way that the corresponding credible regions have
##               exact (or approximately) the nominal frequentist coverage probability.
##
## DEPENDS ON:   'rq' function in package 'quantreg'.
##-----------------------------------------------------------------------------------------------


library(quantreg)


# Simple Metropolis-Hastings

mh <- function(x0, lf, ldprop, rprop, N, B) {

  x <- matrix(NA, N + B, length(x0))
  lfx <- rep(NA, N + B)
  x[1,] <- x0
  lfx[1] <- lf(x0)
  ct <- 0
  for(i in 2:(N + B)) {

    u <- rprop(x[i-1,])
    lfu <- lf(u)
    r <- lfu + ldprop(x[i-1,], u) - lfx[i-1] - ldprop(u, x[i-1,])
    R <- min(exp(r), 1)
    if(ifelse(is.na(R), FALSE, runif(1) <= R)) {

      ct <- ct + 1
      x[i,] <- u
      lfx[i] <- lfu

    } else {

      x[i,] <- x[i-1,]
      lfx[i] <- lfx[i-1]

    }

  }
  x=x[-(1:B),]
  lfx = lfx[-(1:B)]
  return(list(x=x, lfx=lfx))

}


# Main GPS function (for quantile regression, median case)

gps.qr <- function(data, M, B, w, alpha) {

  n <- nrow(data)
  eps <- 1.1 / B
  ldprop <- function(x, x0) sum(dnorm(x, x0, 0.3, log=TRUE))
  rprop <- function(x0) rnorm(2, x0, 0.3)
  theta.hat <- as.numeric(rq(data[,2] ~ data[,1], tau=0.5)$coef)
  data.star <- list()
  cvg <- rep(0, B)
  for(b in 1:B) {

    id <- sample(n, n, replace=TRUE)
    data.star[[b]] <- data[id,]

  }
  t <- 1
  k <- function(t) (1 + t)**(-0.51)
  repeat {

    for(b in 1:B) {

      lf <- function(u) -w * sum(abs(data.star[[b]][,2] - u[1] - u[2] * data.star[[b]][,1])) / 2
      cut <- as.numeric(quantile(mh(theta.hat, lf, ldprop, rprop, M, 100)$lfx, alpha))
      cvg[b] <- as.numeric(lf(theta.hat) >= cut)

    }
    diff <- mean(cvg) - (1 - alpha); print(c(diff + 1 - alpha, w))
    if(abs(diff) <= eps) break else {

      t <- t + 1
      w <- w + k(t) * diff

    }

  }
  lf <- function(u) -w * sum(abs(data[,2] - u[1] - u[2] * data[,1])) / 2
  final.gibbs <- mh(theta.hat, lf, ldprop, rprop, M, 1000)$x
  intv <- get.cr(final.gibbs, alpha)
  return(list(w=w, intv=intv))

}


# Functions for the simulation example in the paper

get.cr.cov <- function(intv, theta) {

  n <- nrow(intv)
  o <- rep(NA, n)
  for(i in 1:n) o[i] <- (intv[i,1] <= theta[i] && intv[i,2] >= theta[i])
  return(all(o))

}


rmodel <- function(n, theta) {

  X <- rchisq(n, 2) - 2
  Y <- theta[1] + theta[2] * X + rnorm(n, 0, 2)
  return(cbind(X, Y))

}


gps.qr.sim <- function(n, reps, M, B, w, alpha=0.05) {

  theta.true <- c(2, 1)
  coverage <- ww <- rep(0, reps)
  len <- matrix(0, reps, 2)
  for(j in 1:reps) {

    cat("Iteration =", j, "\n")
    o <- gps.qr(rmodel(n, theta.true), M, B, w, alpha)
    ww[j] <- o$w
    coverage[j] <- get.cr.cov(o$intv, theta.true)
    len[j,] <- o$intv[,2] - o$intv[,1]

  }
  return(list(w=ww, cvg=mean(coverage), len=apply(len, 2, mean)))

}


## TO RUN:
## system.time(o <- gps.qr.sim(100, 10, 2000, 100, 1)) / 10; print(o)
