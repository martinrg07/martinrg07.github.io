
##-----------------------------------------------------------------------------------------------
## TITLE:        ebreg -- Empirical Bayes model selection and estimation in a sparse high
##               dimensional Gaussian linear regression model.
##
## VERSION:      2nd version (04/10/2018).
##
## AUTHORS:      R. Martin (rgmartin@uic.edu), R. Mess, and S. G. Walker.
##
## REFERENCE:    Martin, Mess, and Walker, "Empirical Bayes posterior concentration in sparse
##               high-dimensional linear models," arXiv:1406.7718
##
## DESCRIPTION:  Method runs Metropolis-Hastings chain to sample from the marginal posterior
##               distribution for the model; if desired, samples of the regression coefficients
##               for the given model can also be obtained.  Variable inclusion probabilities can
##               be directly obtained from the posterior samples of the models which, in turn,
##               can be used to construct a model selection procedure.
##
## DEPENDS ON:   'lars' function in package 'lars'.
##-----------------------------------------------------------------------------------------------


# INPUT:
# y = vector of respons variables
# X = matrix of predictor variables
# sig2 = error variance, if NULL (default), variance is estimated from data
# log.f = log of the prior for the model size; see 'dcomplex' below
# alpha = likelihood fraction
# gam = conditional prior precision parameter
# M = Monte Carlo sample size (burn-in of size 0.2 * M automatically added)
# sample.beta = logical; if TRUE, samples of beta are obtained

# OUTPUT:
# beta = matrix with rows containing sampled beta, if sample.beta=TRUE
# S = matrix with rows containing the sampled models
# s = vector containing the size of the sampled models
# i = vector containing integer labels for the distinct models sampled


ebreg <- function(y, X, sig2=NULL, log.f, alpha, gam, M, sample.beta=FALSE) {

  n <- nrow(X)
  p <- ncol(X)
  if(!exists("lars")) library(lars)
  o.lars <- lars(X, y, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  cv <- cv.lars(X, y, plot.it=FALSE, se=FALSE, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  b.lasso <- coef(o.lars, s=cv$index[which.min(cv$cv)], mode="fraction")
  S <- as.numeric(b.lasso != 0)
  if(is.null(sig2)) {

    z <- as.numeric(y - X[, S > 0] %*% b.lasso[S > 0])
    sig2 <- sum(z**2) / max(n - sum(S), 1)

  }
  v <- 1 + alpha / gam / sig2
  sq.v <- sqrt(v)
  B <- round(0.2 * M)
  bb <- if(sample.beta) matrix(0, nrow=B + M, ncol=p) else NULL
  SS <- matrix(0, nrow=B + M, ncol=p)
  ss <- numeric(B + M)
  ii <- integer(B + M)
  SS[1,] <- S
  ss[1] <- s <- sum(S)
  ii[1] <- i <- iuS <- nuS <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(S, y, X, sample.beta)
  lprior <- log.f(0:n) - lchoose(p, 0:n)
  lpost <- lprior[s] - alpha * out[[1]]$sse / 2 / sig2 - s * log(v) / 2
  for(m in 1:(B + M)) {

    S.new <- rprop(S, n)
    s.new <- sum(S.new)
    i.new <- compare.to.rows(as.matrix(SS[iuS,]), S.new)
    if(i.new == 0) o.new <- get.lm.stuff(S.new, y, X, sample.beta) else o.new <- out[[i.new]]
    lpost.new <- lprior[s.new] - alpha * o.new$sse / 2 / sig2 - s.new * log(v) / 2
    if(runif(1) <= exp(lpost.new - lpost)) {

      S <- S.new
      s <- s.new
      lpost <- lpost.new
      if(i.new == 0) {

        nuS <- nuS + 1
        i <- nuS
        iuS <- c(iuS, m)
        out[[nuS]] <- o.new

      } else i <- i.new

    }
    SS[m,] <- S
    ss[m] <- s
    ii[m] <- i
    if(sample.beta) {

      b <- numeric(p)
      b[S > 0] <- out[[i]]$b.hat + sq.v * out[[i]]$U %*% rnorm(s)
      bb[m,] <- b

    }

  }
  if(sample.beta) bb[-(1:B),]
  return(list(beta=bb, S=SS[-(1:B),], s=ss[-(1:B)], i=ii[-(1:B)]))

}


# Example

dcomplex <- function(x, n, p, a, b, log=TRUE) {

  o <- -x * (log(b) + a * log(p)) + log(x <= n)
  if(!log) o <- exp(o)
  return(o)

}


ebreg.example <- function(n=100, p=500, r=0.25, beta=0.6 * 1:5, M=2000) {

  sample.beta <- FALSE
  sig2 <- 1
  alpha <- 0.999
  gam <- 0.001
  log.f <- function(x) dcomplex(x, n, p, 0.05, 1)
  s0 <- length(beta)
  beta0 <- c(beta, numeric(p - s0))
  S0 <- as.numeric(beta0 != 0)
  sig2.hat <- NULL
  if(r != 0) {

    R <- (1 - r) * diag(p) + array(r, c(p, p))
    e <- eigen(R)
    sqR <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)

  }
  X <- matrix(rnorm(n * p), nrow=n, ncol=p)
  if(r != 0) X <- X %*% sqR
  y <- as.numeric(X[, 1:s0] %*% beta) + sqrt(sig2) * rnorm(n)
  o <- ebreg(y, X, sig2.hat, log.f, alpha, gam, M, sample.beta)
  sz <- sort(unique(o$s))
  pr <- as.numeric(table(o$s)) / M
  incl.pr <- apply(o$S, 2, mean)
  op <- par(pty="m", mfrow=c(2,1), mar=c(4.2, 4.2, 1, 1))
  plot(sz, pr, type="h", xlab="Model Size", ylab="Probability", ylim=c(0,1))
  abline(v=s0, lty=3)
  plot(incl.pr, xlab="Variable Index", ylab="Inclusion Probability", type="h", ylim=c(0,1))
  par(op)
  return(o)

}


# Auxiliary functions called by 'ebreg'

rprop <- function(S, n) {

 s <- sum(S)
 if(s == n) { S[sample(which(S == 1), 1)] <- 0 }
 else if(s  == 1) { S[sample(which(S == 0), 1)] <- 1 }
 else {

   if(runif(1) <= 0.5) S[sample(which(S == 1), 1)] <- 0
   else S[sample(which(S == 0), 1)] <- 1

 }
 return(S)

}


compare.to.rows <- function(SS, S) {

  h <- function(v) sum(abs(v - S))
  o <- apply(SS, 1, h)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])

}


get.lm.stuff <- function(S, y, X, sample.beta) {

  if(!exists("ginv")) library(MASS)
  X.S <- as.matrix(X[, S > 0])
  o <- lm(y ~ X.S - 1)
  sse <- sum(o$residuals**2)
  if(sample.beta) {

    b.hat <- o$coefficients
    XtX.S <- t(X.S) %*% X.S
    V <- ginv(XtX.S)
    U <- chol(V)

  } else {

    U <- NULL
    b.hat <- NULL

  }
  return(list(sse=sse, b.hat=b.hat, U=U))

}




