
## gimcens: Version 1
## R codes associated with the paper "Generalized IMs for censored data"
## Prepared by Joyce Cahoon
## Note: a more user-friendly version is coming soon


# Libraries ---------------------------------------------------------------
library(maxLik)
library(survival)


# Application with data from PBC trial ------------------------------------
data(pbc)
data <- list(
  X = pbc$time[1:312], 
  delta = as.numeric(pbc$status == 0)[1:312])

# MLE estimate
calculate_MLE <- function(X, indicators){
  N <- sum(indicators == 0)
  return(sum(X)/N)
}

data_mle <- calculate_MLE(data$X, data$delta)
mle_ci <- data_mle*exp(c(-1,1)*abs(qnorm(0.05/2))/sqrt(sum(data$delta == 0)))

# IM estimate
thetaTry <- c(
  seq(.9*mle_ci[1], data_mle, len = 50),
  seq(data_mle + 1, 1.1*mle_ci[2], len = 50))

calculate_log_likelihood <- function(X, indicators, theta){
  sum(pexp(X[indicators == 1], rate = 1/theta, lower.tail = FALSE, log.p = TRUE), 
      dexp(X[indicators == 0], rate = 1/theta, log = TRUE))
}

denominator <- calculate_log_likelihood(data$X, data$delta, data_mle)

# Reversed KM estimator 
km_calculator <- function(X, indicators){
  n <- length(X)
  censored <- sort(unique(X[indicators == 1]))
  if(length(censored) > 0){
    weights <- numeric(length(censored))
    for(c in seq_along(censored)){
      matches <- sum(X[indicators == 1] == censored[c])  
      at_risk <- sum(X >= censored[c])
      weights[c] <- 1-matches/at_risk 
    } 
    product_limit <- c(1, cumprod(weights))
    censored_weights <- abs(diff(product_limit))
    if(sum(censored_weights) != 1){
      # add in an upper detection limit 
      censored <- c(censored, Inf)
      censored_weights <- c(censored_weights, 1-sum(censored_weights))
    }
  } else {
    censored <- Inf 
    censored_weights <- 1
  }
  return(list(c = censored, 
              w = censored_weights, 
              p = product_limit))
}

censored <- km_calculator(data$X, data$delta)

# Plausibility values at various theta hypotheses
mc <- 500
n <- length(data$X)
pl <- numeric(length(thetaTry))
for(i in 1:length(thetaTry)){
  theta <- thetaTry[i]
  numerator <- calculate_log_likelihood(data$X, data$delta, theta)
  observed <- exp(numerator-denominator)
  possible_T <- numeric(mc)
  
  for(m in 1:mc){
    new <- rexp(n, rate = 1/theta)
    if(length(censored$c) > 1){
      ci <- sample(censored$c, size = n, replace = TRUE, prob = censored$w)
    } else{
      ci <- rep(censored$c, n)
    }
    new_Y <- pmin(new, ci)
    new_delta <- (new > ci)*1
    new_numerator <- calculate_log_likelihood(new_Y, new_delta, theta)
    new_mle <- calculate_MLE(new_Y, new_delta)
    new_denominator <- calculate_log_likelihood(new_Y, new_delta, new_mle)
    possible_T[m] <- exp(new_numerator-new_denominator)
  }
  pl[i] <- sum(possible_T <= observed)/mc
}  

# Visual
theta_range <- length(pl)
idx <- which.max(pl)
smoothS1 <- smooth.spline(thetaTry[1:idx], pl[1:idx], spar=0.5)
smoothS2 <- smooth.spline(thetaTry[idx:theta_range], pl[idx:theta_range], spar=0.5)
plot(
  thetaTry, pl,
  xlim = c(3500, 5500),
  xlab = expression(theta), 
  ylab = "Possibility", 
  col = "white")
lines(smoothS1, lwd = 2)
lines(smoothS2, lwd = 2)
abline(h = .05, lty = 2, col = "grey5")
segments(mle_ci[1], .025, mle_ci[2], .025, col= 'grey5')
segments(mle_ci[1], .05, mle_ci[1], 0, col= 'grey5')
segments(mle_ci[2], .05, mle_ci[2], 0, col= 'grey5')


# Application with data from ovarian cancer trial -------------------------
rm(list=ls())

data(ovarian)
data <- list(
  X = ovarian$futime, 
  delta = as.numeric(ovarian$fustat == 0))

# MLE estimate 
calculate_MLE <- function(X, indicators, alpha0, beta0, se = FALSE){
  logLikeCensored <- function(params){
    alphahat <- params[1]
    betahat <- params[2]
    if(alphahat < 0 | betahat < 0) return(NA)
    return(sum(pweibull(X[indicators == 1], shape = alphahat, scale = betahat, log.p = TRUE, lower.tail = FALSE),
               dweibull(X[indicators == 0], shape = alphahat, scale = betahat, log = TRUE)))
  }
  fit <- maxLik(logLik = logLikeCensored,
                start = c(alpha0, beta0),
                method = "NM")
  if (se) {
    fit_sd <- stdEr(fit)
    return(list(alphahat = fit$estimate[1],
                betahat = fit$estimate[2], 
                alphahatsd = fit_sd[1], 
                betahatsd = fit_sd[2]))
  } else {
    return(list(alphahat = fit$estimate[1],
                betahat = fit$estimate[2])) 
  }
}

mle <- calculate_MLE(data$X, data$delta, 1, 1, se = TRUE)

# IM estimate
alpha.seq <- seq(.1, 3, len = 50)
beta.seq <- seq(500, 4000, len = 90)
pl <- 0*outer(alpha.seq, beta.seq)
mc <- 500
n <- length(data$X)

calculate_log_likelihood <- function(X, indicators, alpha0, beta0){
  sum(pweibull(X[indicators == 1], shape = alpha0, scale = beta0, log.p = TRUE, lower.tail = FALSE),
      dweibull(X[indicators == 0], shape = alpha0, scale = beta0, log = TRUE))
}

denominator <- calculate_log_likelihood(data$X, data$delta, mle$alphahat, mle$betahat)

km_calculator <- function(X, indicators){
  n <- length(X)
  censored <- sort(unique(X[indicators == 1]))
  if(length(censored) > 0){
    weights <- numeric(length(censored))
    for(c in seq_along(censored)){
      matches <- sum(X[indicators == 1] == censored[c])  
      at_risk <- sum(X >= censored[c])
      weights[c] <- 1-matches/at_risk 
    } 
    product_limit <- c(1, cumprod(weights))
    censored_weights <- abs(diff(product_limit))
    if(sum(censored_weights) != 1){
      # add in an upper detection limit 
      censored <- c(censored, Inf)
      censored_weights <- c(censored_weights, 1-sum(censored_weights))
    }
  } else {
    censored <- Inf # some super low censoring value 
    censored_weights <- 1
  }
  return(list(c = censored, 
              w = censored_weights,
              p = product_limit))
}

cens <- km_calculator(data$X, data$delta)

for(row in 1:50){
  for(col in 1:90){
    alphaTry <- alpha.seq[row]
    betaTry <- beta.seq[col]
    numerator <- calculate_log_likelihood(data$X, data$delta, alphaTry, betaTry)
    logratio <- numerator - denominator
    
    possibleT <- numeric(mc)
    for(m in 1:mc){
      c <- sample(cens$c, size = n, replace = TRUE, prob = cens$w)
      Y <- rweibull(n, shape = alphaTry, scale = betaTry)
      while(sum(Y == 0) > 0){
        Y <- rweibull(n, shape = alphaTry, scale = betaTry)
      }
      X <- pmin(Y, c)
      delta <- (Y > c)*1
      new_numerator <- calculate_log_likelihood(X, delta, alphaTry, betaTry)
      new_mle <- calculate_MLE(X, delta, alphaTry, betaTry)
      new_denominator <- calculate_log_likelihood(X, delta, new_mle$alphahat, new_mle$betahat)
      possibleT[m] <- new_numerator - new_denominator
    }
    pl[row, col] <- sum(possibleT <= logratio) / mc
    cat("Now working on row ", row, " and col ", col, "\n")
  }
}

# Visual 
contour(
  alpha.seq, beta.seq, pl, level = c(.1, .5, .9), 
  xlab = expression(paste("Shape (", beta, ")")),
  ylab = expression(paste("Scale (", lambda, ")")))


# Application in Atrazine concentration data ------------------------------
rm(list=ls())

# Data from Krishnamoorthi et al.
delta <- c(0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1)
X <- c(0.38, 0.05, 0.01, 0.03, 0.03, 0.05, 0.02, 0.01, 0.01, 0.01, 0.11, 0.09, 
       0.01, 0.01, 0.01, 0.01, 0.02, 0.05, 0.02, 0.02, 0.05, 0.03, 0.05, 0.01)

# Transform
Y <- log(X)

# MLE estimate
calculate_MLE <- function(logdata, indicators, mu, sigma2, se = FALSE){
  logLikCensoredFun <- function(params){
    mean <- params[1]
    variance <- params[2]
    if(variance < 0 | is.nan(variance)) return(NA)
    return(sum(pnorm(logdata[indicators == 1], mean, sqrt(variance), log.p = TRUE),
               dnorm(logdata[indicators == 0], mean, sqrt(variance), log = TRUE)))
  }
  fit <- maxLik(logLik = logLikCensoredFun,
                start = c(mu, sigma2),
                method = "NM")
  if (se) {
    fit_sd <- stdEr(fit)
    return(list(muhat = fit$estimate[1], 
                sigma2hat = fit$estimate[2],  
                musdhat = fit_sd[1], 
                sigma2sdhat = fit_sd[2]))
  } else {
    return(list(muhat = fit$estimate[1], 
                sigma2hat = fit$estimate[2]))
  }
}

data_mle <- calculate_MLE(Y, delta, mean(Y), 1, se = TRUE)

# IM estimate
muSeq <- seq(-5, -3, len = 50)
sigSeq <- seq(.2, 6, len= 40)
pl <- 0*outer(muSeq, sigSeq)
mc <- 500
n <- length(Y)  

calculate_log_likelihood <- function(logdata, indicators, mu, sigma2){
  sum(pnorm(logdata[indicators == 1], mean = mu, sd = sqrt(sigma2), log.p = TRUE), 
      dnorm(logdata[indicators == 0], mean = mu, sd = sqrt(sigma2), log = TRUE))
}

denominator <- calculate_log_likelihood(Y, delta, data_mle$muhat, data_mle$sigma2hat)

km_calculator <- function(logdata, indicators){
  n <- length(logdata)
  censored <- sort(unique(logdata[indicators == 1]), decreasing = TRUE)
  if(length(censored) >0){
    weights <- numeric(length(censored))
    for(c in seq_along(censored)){
      matches <- sum(logdata[indicators == 1] == censored[c])  
      at_risk <- sum(logdata <= censored[c])
      weights[c] <- 1-matches/at_risk 
    } 
    product_limit <- c(1, cumprod(weights))
    censored_weights <- abs(diff(product_limit))
    if(sum(censored_weights) != 1){
      # add in an even lower detection limit 
      lowest <- -Inf
      censored <- c(censored, lowest)
      censored_weights <- c(censored_weights, 1-sum(censored_weights))
    }
  } else {
    censored <- -Inf 
    censored_weights <- 1
  }
  return(list(c = censored, 
              w = censored_weights))
}

cens <- km_calculator(Y, delta)

for(row in 1:length(muSeq)){
  for(col in 1:length(sigSeq)){
    muTry <- muSeq[row]
    sTry <- sigSeq[col]
    numerator <- calculate_log_likelihood(Y, delta, muTry, sTry)
    logratio <- numerator - denominator
    
    possibleT <- numeric(mc)
    for(m in 1:mc){
      ci <- sample(cens$c, size = n, replace = TRUE, prob = cens$w)
      new <- rnorm(n, muTry, sqrt(sTry))
      new_Y <- pmax(new, ci)
      new_delta <- (new < ci)*1
      new_numerator <- calculate_log_likelihood(new_Y, new_delta, muTry, sTry)
      new_mle <- calculate_MLE(new_Y, new_delta, muTry, sTry)
      new_denominator <- calculate_log_likelihood(new_Y, new_delta, new_mle$muhat, new_mle$sigma2hat)
      possibleT[m] <- new_numerator - new_denominator
    }
    pl[row, col] <- sum(possibleT <= logratio) / mc
    cat("Now working on row ", row, " and col ", col, "\n")
  }
}

# Visual 
psi <- function(m, s) exp(m + s / 2)
psi.span <- outer(muSeq, sigSeq, psi)
psi.grid <- psi.pl <- c()
pl.levels <- seq(0.001, 1, len=50)
for(j in seq_along(pl.levels)) {
  psi.pl <- c(psi.pl, rep(pl.levels[j], 2))
  psi.grid <- c(psi.grid, range(psi.span[pl >= pl.levels[j]]))
}
ip <- order(psi.grid)
psi.grid <- psi.grid[ip]
psi.pl <- psi.pl[ip]
cut <- max(psi.grid[psi.pl == 1])

plot(
  psi.grid, psi.pl, type="n", 
  ylim=c(0, 1), 
  xlab=expression(psi), 
  ylab="Plausibility", 
  lwd = 3)
lines(psi.grid[psi.grid <= cut], psi.pl[psi.grid <= cut], lwd = 3)
lines(smooth.spline(psi.grid[psi.grid > cut], psi.pl[psi.grid > cut], df=10), lwd = 3)
mle_answer <- exp(data_mle$muhat + data_mle$sigma2hat/2)
abline(h=0.05, lty=3, col="gray30")
abline(v=mle_answer, lty=2, col="red", lwd = 2)

