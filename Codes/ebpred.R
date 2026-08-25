
##-----------------------------------------------------------------------------------------------
## TITLE:        ebpred -- Empirical Bayes prediction in a sparse high dimensional Gaussian 
##               linear regression model using a prior for the error variance.
##
## VERSION:      1st version (02/20/2019).
##
## AUTHORS:      R. Martin (rgmarti3@ncsu.edu), Y. Tang (ytang22@ncsu.edu) 
##
## DESCRIPTION:  Method runs Metropolis-Hastings chain to sample from the marginal posterior
##               distribution for the model and uses an inverse gamma prior for the error 
##               variance for prediction.  
##
## DEPENDS ON:   'lars' function in package 'lars'.
##-----------------------------------------------------------------------------------------------


# INPUT:
# y = vector of response variables
# X = matrix of predictor variables
# XX = vector of new values for the predictor variables to be tested
# alpha = likelihood fraction
# gam = conditional prior precision parameter
# igpar = the parameters for the inverse gamma prior on the error variance
# log.f = log of the prior for the model size; see 'dcomplex' below
# M = Monte Carlo sample size (burn-in of size 0.2 * M automatically added)

# OUTPUT:
# yy = matrix containing predicted reponses
# Yhat = vector containing the predictions for the predictor values tested, XX
# PI.low = low-end of the 95% prediction interval
# PI.high = high-end of the 95% prediction interval

ebpred.prior <- function(y, X, XX, alpha, gam, igpar, log.f, M) {

  n <- nrow(X)
  p <- ncol(X)
  d <- nrow(XX)
  
  if(!exists("lars")) library(lars)
  o.lars <- lars(X, y, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  cv <- cv.lars(X, y, plot.it=FALSE, se=FALSE, normalize=FALSE, intercept=FALSE, use.Gram=FALSE)
  b.lasso <- coef(o.lars, s=cv$index[which.min(cv$cv)], mode="fraction")
  
  S <- as.numeric(b.lasso != 0)
  v <- 1 + alpha / gam
  B <- round(0.2 * M)
  SS <- matrix(0, nrow=B + M, ncol=p)
  ss <- numeric(B + M)
  ii <- integer(B + M)
  YY <- prediction <- matrix(0, nrow=M, ncol=d)
  SS[1,] <- S
  ss[1] <- s <- sum(S)
  ii[1] <- i <- iuS <- nuS <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(S, y, X)
  lprior <- log.f(0:n) - lchoose(p, 0:n)
  lpost <- lprior[s+1] - s * log(v) / 2 - (alpha * n / 2 + igpar[1]) * log(igpar[2] + alpha * out[[1]]$sse / 2)
  for(m in 1:(B + M)) {

    S.new <- rprop(S, n)
    s.new <- sum(S.new)
    i.new <- compare.to.rows(as.matrix(SS[iuS,]), S.new)
    if(i.new == 0) o.new <- get.lm.stuff(S.new, y, X) else o.new <- out[[i.new]]
    lpost.new <- lprior[s.new+1] - s.new * log(v) / 2 - (alpha * n / 2 + igpar[1]) * log(igpar[2] + alpha * o.new$sse / 2)
    r <- exp(lpost.new-lpost)
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
    if(m > B) {

      g.star <- igpar[1] + alpha * n / 2
      h.star <- igpar[2] + alpha * out[[i]]$sse / 2
      if(sum(S)==0 | is.null(out[[i]]$b.hat)) YY[m - B,] <- sqrt(h.star / g.star) * rt(d, df=2 * g.star)  else {

        if(length(out[[i]]$b.hat)==1) out.mean <- XX[,S>0] * out[[i]]$b.hat else out.mean <- XX[,S > 0] %*% out[[i]]$b.hat
        prediction[m - B,] <- out.mean
        if(d == 1) mat <- matrix(XX[,S > 0], nrow=1) else mat <- XX[,S > 0]
        out.scale <- diag(d) + mat %*% out[[i]]$U %*% t(out[[i]]$U) %*% t(mat) / (alpha + gam)
        YY[m - B,] <- out.mean + sqrt(h.star / g.star) * chol(out.scale) %*% rt(d, df=2 * g.star)

      }
    }
  }
  PI <- apply(YY, 2, function(x) quantile(x, probs=c(.025,.975)))
  return(list(yy=YY, Yhat=prediction, PI.low=PI[1], PI.high=PI[2]))

}


# Example...

dcomplex <- function(x, n, p, a, b) -x * (log(b) + a * log(p)) + log(x <= n)

ebpred.sim <- function(reps=100, n=70, p=100, signal=1, r=0.5, M=5000, alpha=.99, gam=.005) {

  beta <- c(0,0,rep(signal,2),rep(0,10),signal,rep(0,6),signal,0,0,signal)
  beta.index <- c(3,4,15,22,25)
  sig2 <- 1
  d <- 1
  log.f <- function(x) dcomplex(x, n, p, 0.05, 1)
  s0 <- length(beta.index)
  g <- function(i, j) r**(abs(i - j))
  R <- outer(1:p, 1:p, g)
  e <- eigen(R)
  sqR <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  time.eb <- mspe.eb <- PIhigh.eb <- PIlow.eb <- trueY.eb <- oracle.low <- oracle.high <- numeric(reps)

  for(k in 1:reps)  {

    X <- matrix(rnorm(n * p), nrow=n, ncol=p) %*% sqR
    X.new <- matrix(rnorm(p), nrow=1, ncol=p) %*% sqR
    y <- as.numeric(X[, beta.index] %*% beta[beta.index]) + sqrt(sig2) * rnorm(n)
    y.new <- as.numeric(X.new[,beta.index] %*% beta[beta.index]) + sqrt(sig2) * rnorm(1)
    out <- lm(y ~ X[,beta.index] - 1)
    or.se <- summary(out)$sigma * sqrt(1 + t(X.new[beta.index]) %*% solve(t(X[,beta.index])%*%X[,beta.index]) 
                                       %*% X.new[beta.index])
    yhat <- t(X.new[beta.index]) %*% out$coefficients
    stddev <- qt(0.975, df=n - 5) * or.se

    oracle.high[k] <- yhat + stddev
    oracle.low[k] <- yhat - stddev
    time.eb[k] <- system.time(o <- ebpred.prior(y, X, X.new, alpha, gam, c(0.01, 4), log.f, M))[3]

    Y.hat <- colMeans(o$Yhat)

    if(reps == 1) {

      op <- par(pty="m", mar=c(4.2, 4.2, 1, 1))
      hist(o$yy, freq=FALSE, breaks=25, xlab=expression(tilde(y)), col="gray", border="white", main="")
      or.f <- function(z) dt((z - Y.hat) / c(or.se), df=n-5) / c(or.se)
      curve(or.f, add=TRUE)
      par(op)

    }
    
    mspe.eb[k] <- mean((y.new - Y.hat)**2) #mspe for empirical Bayes
    PIhigh.eb[k] <- o$PI.high #prediction interval high-end for empirical Bayes
    PIlow.eb[k] <- o$PI.low #prediction interval low-end for empirical Bayes

    if(PIlow.eb[k]<y.new && y.new<PIhigh.eb[k]) trueY.eb[k] <- 1 #coverage of PI for empirical Bayes

  }
  return(cbind(n, p, r, signal, alpha, gam, reps, mspe.eb=mean(mspe.eb), time.eb=mean(time.eb), sd.eb=sd(mspe.eb), 
               len.eb=mean(PIhigh.eb - PIlow.eb), cvg.eb=sum(trueY.eb)/reps, len.or=mean(oracle.high - oracle.low)))
}



# Auxiliary functions called by 'ebpred.prior'


rprop <- function(S, n) {
  
  s <- sum(S)
  if(s == n) { S[sample(which(S == 1), 1)] <- 0 }
  else if(s  == 0) { S[sample(which(S == 0), 1)] <- 1 }
  else if(s  == 1) { S[sample(which(S == 0), 1)] <- 1 }
  else {
    
    if(runif(1) <= 0.5) S[sample(which(S == 1), 1)] <- 0
    else S[sample(which(S == 0), 1)] <- 1
    
  }
  return(S)
  
}


compare.to.rows <- function(SS, S) {
  if (ncol(SS)==1){
    SS <- t(SS)
  }
  numrows <- nrow(SS)
  newS <- matrix(rep(S,numrows), nrow=numrows, byrow=TRUE)
  subtracted <- abs(newS-SS)
  o <- rowSums(subtracted)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])
}


get.lm.stuff <- function(S, y, X) {
  
  if(sum(S)==0){
    return(list(sse=sum(y**2), b.hat=NULL, U=NULL))
  }
  else{
    X.S <- as.matrix(X[, S > 0])
    o <- lm.fit(X.S, y, singular.ok = FALSE)
    sse <- sum(o$residuals**2)
    b.hat <- o$coefficients
    V <- chol2inv(qr.R(o$qr))
    U <- chol(V)
    return(list(sse=sse, b.hat=b.hat, U=U))
  }
  
}

get.mode <- function(x) {
  
  # x is an integer vector
  xt <- table(x)
  mode <- as.integer(names(sort(-xt)[1]))
  return(mode)
  
}
