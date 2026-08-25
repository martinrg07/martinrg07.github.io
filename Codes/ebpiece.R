
##-----------------------------------------------------------------------------------------------
## TITLE:        ebpiece -- Empirical Bayes inference on a piecewise constant normal mean
##
## VERSION:      1st version (12/12/2017).
##
## AUTHORS:      R. Martin (rgmarti3@ncsu.edu) and W. Shen.
##
## REFERENCE:    Martin and Shen, "Asymptotically optimal empirical Bayes inference in a piecewise 
##               constant sequence model," arXiv:...
##
## DESCRIPTION:  Method runs Metropolis-Hastings chain to sample from the marginal posterior
##               distribution for the block configuration; given the block configuration, it is 
##               easy to sample from the normal posterior for in-block means.
##-----------------------------------------------------------------------------------------------


# INPUT:
# y = vector of observations
# sig2 = variance
# alpha = likelihood fraction
# v = prior variance
# lambda = block config prior parameter
# M = Monte Carlo sample size (burn-in of size 0.5 * M automatically added)

# OUTPUT:
# theta = matrix with rows containing sampled theta n-vectors
# B = matrix with rows containing the sampled block configurations
# b = vector containing the size of the sampled block configurations

ebpiece <- function(y, sig2=1, alpha, v, lambda, M) {
  
  n <- length(y)
  burn <- round(0.5 * M)
  MM <- M + burn
  BB <- matrix(0, nrow=MM, ncol=n - 1)
  theta <- matrix(0, nrow=MM, ncol=n)
  B <- sample(0:1, n - 1, replace=TRUE)
  b <- sum(B) + 1
  bb <- numeric(MM)
  BB[1,] <- B
  bb[1] <- b
  ii <- integer(MM)
  ii[1] <- i <- iuB <- nuB <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(B, y)
  vv <- 1 + alpha * v / sig2 
  sq.vv <- sqrt(vv)
  lprior <- -lambda * (1:n - 1) * log(n) - lchoose(n - 1, 1:n - 1)
  lpost <- lprior[b] - alpha * out[[1]]$sse / 2 / sig2 - b * log(vv) / 2
  for(m in 1:MM) {
    
    B.new <- rprop(B)
    b.new <- sum(B.new) + 1
    i.new <- compare.blocks(as.matrix(BB[iuB,]), B.new)
    if(i.new == 0) o.new <- get.lm.stuff(B.new, y) else o.new <- out[[i.new]]
    lpost.new <- lprior[b.new] - alpha * o.new$sse / 2 / sig2 - b.new * log(vv) / 2
    if(runif(1) <= exp(lpost.new - lpost)) {
      
      B <- B.new
      b <- b.new
      lpost <- lpost.new
      if(i.new == 0) {
        
        nuB <- nuB + 1
        i <- nuB
        iuB <- c(iuB, m)
        out[[nuB]] <- o.new
        
      } else i <- i.new
      
    }
    BB[m,] <- B
    bb[m] <- b
    ii[m] <- i
    theta[m,] <- rep(out[[i]]$theta.B + sqrt(v / vv / out[[i]]$B.split) * rnorm(b), out[[i]]$B.split)
    
  }
  return(list(theta=theta[-(1:burn),], B=BB[-(1:burn),], b=bb[-(1:burn)]))
  
  
}


# Auxiliary functions

rprop <- function(B) {
  
  i <- sample(seq_along(B), 1)
  B[i] <- (B[i] + 1) %% 2
  return(B)
  
}


compare.blocks <- function(B.old, B.new) {
  
  h <- function(v) sum(abs(v - B.new))
  o <- apply(B.old, 1, h)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])
  
}


get.lm.stuff <- function(B, y) {
  
  if(all(B == 0)) {
    
    theta.B <- mean(y)
    sse <- sum((y - theta.B)**2)
    B.split <- length(y)
    
  }
  else {
  
    x <- as.factor(c(1, 1 + cumsum(B)))
    o <- lm(y ~ x - 1)
    sse <- sum(o$residuals**2)
    theta.B <- o$coefficients
    B.split <- as.numeric(table(x))
    
  }
  return(list(sse=sse, theta.B=theta.B, B.split=B.split))
  
}


# Do examples

ebpiece.example <- function(theta, sig2, v, lambda, M) {

  n <- length(theta)
  y <- theta + sqrt(sig2) * rnorm(n)
  o <- ebpiece(y, sig2, 0.99, v, lambda, M)
  plot(y, cex=0.5, col="gray")
  lines(theta)
  points(apply(o$theta, 2, mean), pch=19, cex=0.5)
  tb <- table(o$b)
  plot(tb / sum(tb), xlab="|B|", ylab="Probability")
  return(o)

}


if(FALSE) {
  
  
  # Example 1 in paper, from Frick, Munk, and Seiling (2014), p 516
  
  n <- 497
  theta <- numeric(n)
  theta[1:138] <- -0.18
  theta[139:225] <- 0.08
  theta[226:242] <- 1.07
  theta[243:299] <- -0.53
  theta[300:308] <- 0.16
  theta[309:332] <- -0.69
  theta[333:n] <- -0.16
  o <- ebpiece.example(theta, 0.04, 1, 1, 10000)
  
  # Example 2 in paper
  
  n <- 1000
  b0 <- 20
  mag <- 2
  utheta <- runif(b0, -mag, mag)
  theta <- rep(utheta, rep(n / b0, b0)); plot(theta, type="l")
  o <- ebpiece.example(theta, 0.5**2, 1, 1, 10000)
  
}



