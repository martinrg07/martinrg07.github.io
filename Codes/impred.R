
##-----------------------------------------------------------------------------------
## TITLE:      R codes for inferential model (IM) prediction examples
##
## AUTHORS:    R. Martin (rgmartin@uic.edu) and R. Lingham
##
## REFERENCE:  Martin and Lingham, "Prior-free probabilistic prediction of future
##             observations", arXiv:
##-----------------------------------------------------------------------------------


## Log-normal example

lognorm.pred <- function(Y, M, nmean) {

  X <- log(Y)
  Xbar <- mean(X)
  S <- sd(X)
  n <- length(Y)
  f <- function(u) {

    u1 <- u[1]
    u2 <- u[2]
    uu <- u[3:(2 + nmean)]
    d <- (uu - u1 / sqrt(n)) / u2
    return(mean(exp(Xbar - S * d)))

  }
  U <- matrix(rnorm(M * (nmean + 2)), nrow=M, ncol=nmean + 2)
  U[,2] <- sqrt(rchisq(M, df=n - 1) / (n - 1))
  YY <- apply(U, 1, f)
  return(YY)

}


lognorm.pred.sim <- function(reps, n, theta, M, nmean=1) {

  set.seed(7)
  mu <- theta[1]
  sig <- sqrt(theta[2])
  o <- numeric(reps)
  for(r in 1:reps) {

    Y <- exp(rnorm(n, mu, sig))
    YY <- lognorm.pred(Y, M, nmean)
    o[r] <- mean(YY <= mean(exp(rnorm(nmean, mu, sig))))

  }
  return(sort(o))

}


# Creates plot in Figure 1(a)

bhaumik.example <- function(yy, M) {

  nmean <- 5
  alpha <- 0.05
  Y <- c(23.0, 63.0, 3.0, 70.0, 16.0, 5.0, 1.0, 57.0, 5.0, 3.0, 24.0, 2.0, 1.0, 48.0, 3.0)
  YY <- lognorm.pred(Y, M, nmean)
  D <- outer(YY, yy, "-") >= 0
  pl <- apply(D, 2, mean)
  plot(yy, pl, type="l", ylim=c(0,1), xlab=expression(tilde(y)), ylab=expression(pl[y](tilde(y))))
  yupper <- as.numeric(quantile(YY, 1 - alpha))
  abline(h=alpha, v=yupper, lty=3)
  print(yupper)
  return(YY)

}


# Creates plot in Figure 1(b)

bhaumik.example.sim <- function() {

  n <- 15
  nmean <- 5
  mu <- 2.167
  sig2 <- 2.3808
  reps <- 5000
  M <- 50000
  o.bhaumik <- lognorm.pred.sim(reps, n, c(mu, sig2), M, nmean)
  x <- (1:reps - 1) / reps
  plot(x, o.bhaumik, type="l", col="gray", xlab=expression(G[Y](tilde(Y))), ylab="CDF")
  abline(a=0, b=1)

}



## Gamma example -- Section 4.3

# Solving for gamma parameters; see Appendix

gamma2.find.par <- function(u, y) {

  n <- length(y)
  t1 <- sum(y)
  t2 <- log(t1 / n) - mean(log(y))
  u1 <- u[1]
  u2 <- u[2]
  g <- function(w) {

    m <- log(w) - digamma(w + 1) + 1 / w              # mean approximation
    v <- (trigamma(w + 1) + 1 / w**2 - 1 / w) / n     # variance approximation
    return(pgamma(t2, m**2 / v, m / v) - u2)          # moment-matching gamma cdf

  }
  a <- 5e-3 / n
  b <- n
  oops <- FALSE
  if(g(a) >= 0) { theta1 <- a; oops <- TRUE }
  b <- n
  if(g(b) <= 0 && !oops) {

    repeat {

      b <- 10 * b
      if(g(b) > 0) break

    }

  }
  if(!oops) theta1 <- uniroot(g, interval=c(a, b))$root
  theta2 <- t1 / qgamma(u1, n * theta1)
  return(c(theta1, theta2))

}


gamma2.pred <- function(Y, M, nmax) {

  rgamma.max <- function(tt) {

    if(tt[1] > 0.001) o <- rgamma(nmax, shape=tt[1]) else o <- exp(-rexp(nmax) / tt[1])
    return(tt[2] * max(o))

  }
  U <- matrix(runif(2 * M), ncol=2)
  TT <- t(apply(U, 1, gamma2.find.par.c, y=Y, UU=runif(5000 * length(Y))))
  YY <- apply(TT, 1, rgamma.max)
  return(YY)

}


gamma2.pred.sim <- function(reps, n, theta, M, nmax) {

  set.seed(7)
  theta1 <- theta[1]
  theta2 <- theta[2]
  o <- numeric(reps)
  for(r in 1:reps) {

    Y <- theta2 * rgamma(n, theta1)
    YY <- gamma2.pred(Y, M, nmax)
    o[r] <- mean(YY <= theta2 * max(rgamma(nmax, shape=theta1)))

  }
  return(o)

}


# Creates plot in Figure 2(a)

hamada.example <- function(yy=50:200, M) {

  Y <- c(18, 23, 29, 409, 24, 74, 13, 62, 46, 4, 57, 19, 47, 13, 19, 208, 119, 209, 10, 188)
  nmax <- 5
  set.seed(77)
  YY <- gamma2.pred(Y, M, nmax)
  D <- outer(YY, yy, "-") <= 0
  pl <- apply(D, 2, mean)
  plot(yy, pl, type="l", ylim=c(0,1), xlab=expression(tilde(y)), ylab=expression(pl[y](tilde(y))))
  ylower <- as.numeric(quantile(YY, 0.1))
  abline(h=0.1, v=ylower, lty=3)
  print(ylower)
  return(list(yy=yy, pl=pl, YY=sort(YY)))

}


# Creates plot in Figure 2(b)

hamada.example.sim <- function(reps=2000, M=10000) {

  n <- 20
  theta <- c(0.8763, 1 / 0.0110)
  nmax <- 5
  o <- gamma2.pred.sim(reps, n, theta, M, nmax)
  x <- (1:reps - 1) / reps
  plot(x, sort(o), type="l", col="gray", xlab=expression(G[Y](tilde(Y))), ylab="CDF")
  abline(a=0, b=1)
  return(o)

}


## Binomial example -- Section 4.4

binom.pred <- function(Y, n, m, M) {

  f <- function(u) {

    if(Y == 0) theta1 <- 0 else theta1 <- qbeta(u[1], Y, n - Y + 1)
    if(Y == n) theta2 <- 1 else theta2 <- qbeta(u[1], Y + 1, n - Y)
    prob <- rep(runif(1, theta1, theta2), 2)
    return(qbinom(u[2], size=m, prob=prob))

  }
  U <- matrix(runif(2 * M), ncol=2)
  YY <- t(apply(U, 1, f))
  return(YY)

}


binom.pred.fi <- function(Y, n, m, M) {

  if(Y == 0) {

    lower <- numeric(M)
    upper <- rbeta(M, 1, n)

  } else if(Y == n) {

    upper <- 1 + numeric(M)
    lower <- rbeta(M, n, 1)

  } else {

    upper <- rbeta(M, Y + 1, n - Y)
    lower <- upper * rbeta(M, Y, 1)

  }
  TT <- runif(M, lower, upper)
  YY <- rbinom(M, size=m, prob=TT)
  return(YY)

}


binom.pred.ba <- function(Y, n, m, M) {

  ay <- 0.5 + Y
  by <- 0.5 + n - Y
  th <- rbeta(M, ay, by)
  YY <- rbinom(M, size=m, prob=th)
  return(YY)

}


# Simulation for plots in Figure 3

binom.pred.sim <- function(reps, n, theta, M, m) {

  alpha <- 0.05
  cvg.im <- cvg.fi <- cvg.ba <- len.im <- len.fi <- len.ba <- 0 * theta
  for(i in seq_along(theta)) {

    set.seed(7)
    p <- theta[i]
    Y <- rbinom(reps, size=n, prob=p)
    Ytil <- rbinom(reps, size=m, prob=p)
    for(r in 1:reps) {

      YY.im <- binom.pred(Y[r], n, m, M)[,2]
      YY.fi <- binom.pred.fi(Y[r], n, m, M)
      YY.ba <- binom.pred.ba(Y[r], n, m, M)
      o.im <- as.numeric(quantile(YY.im, 1-alpha))
      o.fi <- as.numeric(quantile(YY.fi, 1-alpha))
      o.ba <- as.numeric(quantile(YY.ba, 1-alpha))
      len.im[i] <- len.im[i] + o.im / reps
      len.fi[i] <- len.fi[i] + o.fi / reps
      len.ba[i] <- len.ba[i] + o.ba / reps
      cvg.im[i] <- cvg.im[i] + (o.im >= Ytil[r]) / reps
      cvg.fi[i] <- cvg.fi[i] + (o.fi >= Ytil[r]) / reps
      cvg.ba[i] <- cvg.ba[i] + (o.ba >= Ytil[r]) / reps

    }

  }
  o <- cbind(theta, cvg.im, cvg.fi, cvg.ba, len.im, len.fi, len.ba)
  return(o)

}



