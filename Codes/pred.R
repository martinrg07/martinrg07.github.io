
## R code for the simulation experiments carried out in Section 4.2 of 
##
##   Hong, Kuffner, and Martin (201x). "On prediction of future insurance claims
##   when the model is uncertain," submitted for publication; SSRN Abstract 2883574
##
## Source the file and run the code inside if(FALSE) {...} below.  Results summarize 
## the coverage and length of prediction intervals based on the oracle/true, full, 
## AIC-selected models, along with that based on a bootstrap adjustment to the naive 
## plug-in AIC-based intervals.



library(leaps)
library(MASS)



# Simulate data from model Y = X beta + epsilon, where X ~ N(0, R), R is a 
# correlation matrix, and epsilon ~ N(0, sigma^2)

simu <- function(n, beta, sigma, corr.X) {
  
  p <- dim(beta)[1];
  H <- t(chol(corr.X))
  X <- matrix(rnorm(n*p), ncol=p) %*% t(H)
  Y <- X %*% beta + sigma * rnorm(n)
  return(list(X=X, Y=Y))
  
}


# Get predictive distribution and prediction interval based on OLS

ols.pred <- function(X, Y, x, alpha, intercept=FALSE) {
  
  n <- length(Y)
  p <- intercept + ncol(X)
  if(intercept) XX <- cbind(1 + 0 * Y, X) else XX <- X
  M <- ginv(t(XX) %*% XX) 
  beta.hat <- M %*% t(XX) %*% Y
  sigma.hat <- sqrt(sum((Y - XX %*% beta.hat)**2) / (n - p))
  if(intercept) xx <- c(1, x) else xx <- c(x)
  m <- t(xx) %*% beta.hat
  pe <- sigma.hat * sqrt(1 + t(xx) %*% M %*% xx)
  dpred <- function(y) dt((y - m) / pe, df=n - p) / pe
  int <- m + c(-1, 1) * qt(1 - alpha / 2, df=n - p) * pe
  return(list(m=m, sigma.hat=sigma.hat, dpred=dpred, int=int))
  
}


# Get AIC selected variables, etc

aic.select <- function(X, Y, x, intercept) {
  
  n <- length(Y)
  o <- regsubsets(x=X, y=Y, intercept=intercept, nbest=1, nvmax=ncol(X), method="exhaustive")
  o.sum <- summary(o)
  S <- o.sum$which[which.min(o.sum$cp), ]
  if(intercept) {
    
    XX <- cbind(rep(1, n), X)[, S]
    xx <- c(1, x)[S]
    
  } else {
    
    XX <- X[, S]
    xx <- c(x)[S]
    
  }
  return(list(S=S, XX=XX, xx=xx))
  
}


# Get predictive distribution and prediction interval based on AIC
# (These are based on naive application of OLS theory after selection)

aic.pred <- function(X, Y, x, alpha, intercept=FALSE) {
  
  n <- length(Y)
  A <- aic.select(X, Y, x, intercept)
  o.lm <- lm(Y ~ A$XX - 1)
  beta.hat <- o.lm$coefficients
  sig.hat <- summary(o.lm)$sigma
  df <- n - sum(A$S)
  y.hat <- t(A$xx) %*% beta.hat 
  pe <- c(sig.hat * sqrt(1 + t(A$xx) %*% solve(t(A$XX) %*% A$XX) %*% A$xx))
  pred.int <- c(y.hat - pe * qt(1 - alpha / 2, df=df), y.hat + pe * qt(1 - alpha / 2, df=df))
  dpred <- function(y) dt((y - y.hat) / pe, df=df) / pe
  return(list(int=pred.int, S=A$S, y.hat=y.hat, sigma.hat=sig.hat, pe=pe, dpred=dpred))
  
}


# Get predictive distribution and prediction interval based on bootstrap adjusted AIC

aic.boot.pred <- function(X, Y, x, B, alpha, intercept=FALSE) {
  
  n <- length(Y)
  A <- aic.select(X, Y, x, intercept)
  oo <- lm(Y ~ A$XX - 1)
  Yfit <- fitted.values(oo)
  Yhat <- c(t(A$xx) %*% oo$coefficients)
  e <- summary(oo)$sigma * stdres(oo)
  me <- mean(e)
  y.hat <- pm <- ps <- numeric(B)
  for(b in 1:B) {
    
    ee <- sample(e - me, size=n, replace=TRUE)
    YY <- Yfit + ee
    Ab <- aic.select(X, YY, x, intercept)
    ob <- lm(YY ~ Ab$XX - 1)
    pm[b] <- t(Ab$xx) %*% ob$coefficients
    ps[b] <- summary(ob)$sigma 
    e.star <- ps[b] * stdres(ob)
    y.hat[b] <-  2* Yhat - pm[b] + sample(e.star - mean(e.star), size=1)
    
  }
  pred.int <- as.numeric(quantile(y.hat, c(alpha / 2, 1 - alpha / 2)))
  return(list(int=pred.int, y.hat=y.hat, pm=pm, ps=ps))
  
}



# Simulation to assess prediction interval quality 

aic.sim <- function(n, reps, nbeta, p, sigma, rho, rcoef, alpha, B, intercept) {
  
  corr.X <- rho**(abs(outer(1:p, 1:p, "-")))
  cvg.aic <- cvg.aic.b <- cvg.full <- cvg.or <- numeric(nbeta)
  len.aic <- len.aic.b <- len.full <- len.or <- numeric(nbeta)
  for(j in 1:nbeta) {
      
    stuff <- rcoef(p)
    beta <- stuff$beta
    S <- stuff$S
    for(r in 1:reps) {
        
      o <- simu(n, beta, sigma, corr.X)
      o.new <- simu(1, beta, sigma, corr.X)
      aic <- aic.pred(o$X, o$Y, o.new$X, alpha, intercept)
      aic.b <- aic.boot.pred(o$X, o$Y, o.new$X, B, alpha, intercept)
      full <- ols.pred(o$X, o$Y, c(o.new$X), alpha, intercept)
      or <- ols.pred(o$X[, S], o$Y, c(o.new$X[, S]), alpha, intercept)
      aic.int <- aic$int
      aic.b.int <- aic.b$int
      full.int <- full$int
      or.int <- or$int
      cvg.aic[j] <- cvg.aic[j] + (aic.int[1] <= o.new$Y && aic.int[2] >= o.new$Y) / reps
      cvg.aic.b[j] <- cvg.aic.b[j] + (aic.b.int[1] <= o.new$Y && aic.b.int[2] >= o.new$Y) / reps
      cvg.full[j] <- cvg.full[j] + (full.int[1] <= o.new$Y && full.int[2] >= o.new$Y) / reps
      cvg.or[j] <- cvg.or[j] + (or.int[1] <= o.new$Y && or.int[2] >= o.new$Y) / reps
      len.aic[j] <- len.aic[j] + (aic.int[2] - aic.int[1]) / reps
      len.aic.b[j] <- len.aic.b[j] + (aic.b.int[2] - aic.b.int[1]) / reps
      len.full[j] <- len.full[j] + (full.int[2] - full.int[1]) / reps
      len.or[j] <- len.or[j] + (or.int[2] - or.int[1]) / reps
        
    }
      
  }
  cvg <- cbind(aic=cvg.aic, boot=cvg.aic.b, full=cvg.full, or=cvg.or)
  len <- cbind(aic=len.aic, boot=len.aic.b, full=len.full, or=len.or)
  return(list(cvg=cvg, len=len))
  
}


if(FALSE) {
  
  # Run the AIC simulation
  
  n <- 50
  reps <- 200 
  nbeta <- 20
  p <- 10
  sigma <- 1
  rho <- 0.5
  alpha <- 0.05
  B <- 500
  intercept <- TRUE
  signal.size <- 1  #3  #5
  rcoef <- function(P) {
    
    beta <- numeric(P)
    beta[sample(P, 3, replace=FALSE)] <- 1
    S <- as.logical(beta)
    beta[S] <- rnorm(3, signal.size, 1)
    return(list(S=S, beta=matrix(beta, ncol=1)))
    
  }
  out <- aic.sim(n, reps, nbeta, p, sigma, rho, rcoef, alpha, B, intercept)
  boxplot(out$cvg, ylab="Coverage", col="gray", main="")
  abline(h=1-alpha, lty=3)
  boxplot(out$len[,-4] / out$len[,4], ylab="Relative-to-oracle length", col="gray", main="")
  abline(h=1, lty=3)
  
}



