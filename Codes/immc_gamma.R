
## R code for the stitched IM : gamma application
## Ryan Martin (rgmarti3@ncsu.edu)
## 06 Feb 2025
## Reference: https://arxiv.org/abs/2501.10585


# Auxiliary function : evaluate xi(data, alpha)

imvar <- function(xi, alpha, pl, mle, J, tol=1e-2, a=1, b=1, max.it=25) {
  
  xi <- log(xi) 
  D <- length(mle)
  if(D == 1) posts <- sqrt(qchisq(1 - alpha, 1) / J)
  else posts <- t(J$vectors) * sqrt(qchisq(1 - alpha, D) / J$values)
  maxpl <- function(v) max( c(pl(mle + v), pl(mle - v)) )
  w <- function(s) a / (1 + s)**b
  it <- 1
  s.time <- Sys.time()
  repeat {
    
    posts.xi <- posts * exp(xi / 2)
    if(D == 1) g.xi <- maxpl(posts.xi) - alpha
    else g.xi <- apply(posts.xi, 1, maxpl) - alpha
    if(all(abs(g.xi) <= tol) || (it >= max.it)) break else {
      
      xi <- xi + w(it) * g.xi
      it <- it + 1
      
    }
    
  }
  e.time <- Sys.time()
  timedf <- as.numeric(e.time - s.time)
  if(D == 1) iJ.xi <- sqrt(exp(xi) / J)
  else iJ.xi <- J$vectors %*% diag(sqrt(exp(xi) / J$values)) 
  return(list(xi=exp(xi), g.xi=g.xi, iter=it, time=timedf, sqiJ.xi=iJ.xi))
  
}

# Libraries for running parallel computations

library(foreach)
library(parallel)
library(doParallel)


# Specific functions related to the gamma model 

gamma.mle <- function(X) {
  
  n <- length(X)
  T1 <- sum(log(X))
  T2 <- sum(X)
  mme <- c(mean(X)**2 / var(X), var(X) / mean(X))
  #f <- function(th) -sum(dgamma(X, shape=th[1], scale=th[2], log=TRUE))
  lf <- function(lth) -sum(dgamma(X, shape=exp(lth[1]), scale=exp(lth[2]), log=TRUE))
  df <- function(th) {
    
    alpha <- th[1]
    beta <- th[2]
    o1 <- -n * digamma(alpha) - n * log(beta) + T1
    o2 <- -n * alpha / beta + T2 / beta**2
    return(-c(o1, o2))
    
  }
  ddf <- function(th) {
    
    alpha <- th[1]
    beta <- th[2]
    o11 <- -n * trigamma(alpha)
    o22 <- n * alpha / beta**2 - 2 * T2 / beta**3
    o12 <- -n / beta
    return(matrix(-c(o11, o12, o12, o22), nrow=2, byrow=TRUE))
    
  }
  #o <- optim(par=mme, fn=f, gr=df, method="BFGS") 
  o <- optim(par=log(mme), fn=lf, method="BFGS") 
  #mle <- o$par 
  mle <- exp(o$par)
  val <- -o$value
  V <- diag(mle) %*% ddf(mle) %*% diag(mle)
  return(list(mle=log(mle), val=val, obsinf=V))
  
}


logrellik <- function(th, x) {
  
  v <- exp(th)
  lnum <- sum(dgamma(x, shape=v[1], scale=v[2], log=TRUE))
  lden <- gamma.mle(x)$val
  return(lnum - lden)
  
}


gamma.pl <- function(x, theta1, theta2, UU) {
  
  f <- function(y) logrellik(c(theta1, theta2), y)
  f.x <- f(x)
  f.X <- apply(exp(theta2) * qgamma(UU, shape=exp(theta1)), 1, f)
  pl <- mean(f.X <= f.x)
  return(pl)
  
}


# Data and summaries

X <- c(152, 152, 115, 109, 137, 88, 94, 77, 160, 165, 125, 40, 128, 123, 136, 101, 62, 153, 83, 69)
n <- length(X)
  
mle.out <- gamma.mle(X)
mle <- mle.out$mle
J <- mle.out$obsinf
eJ <- eigen(J)
  

# evaluate the possibilistic IM contour 

M <- 500
UU <- matrix(runif(n * M), nrow=M)
pl <- function(z) gamma.pl(X, z[1], z[2], UU)


# evaluate xi(data, alpha) on a grid of alpha's

B <- 100
AA <- seq(0.001, 0.999, length=B)
h <- function(s) {
    
  o <- imvar(1, AA[s], pl, mle, eJ, tol=5e-2, a=5, b=1, max.it=25)
  return(o$xi)
    
}
start.time <- Sys.time()
numCores <- detectCores() 
clus <- makeCluster(numCores)
registerDoParallel(clus)
Xi <- t( foreach(b=1:B, .combine='cbind') %dopar% { h(b) } )
stopCluster(clus)
end.time <- Sys.time()
print(as.numeric(end.time - start.time))
  

# Get samples from the Q^* approximation

N <- 5000
A <- runif(N) 
Z <- matrix(rnorm(2 * N), nrow=N)
H <- function(k) {
    
  a <- A[k]
  u <- Z[k,] / sqrt(sum(Z[k,]**2))
  if(any(AA == a)) xi.a <- Xi[which(AA == a),] 
  else if(a < min(AA)) xi.a <- Xi[1,] #jitter(Xi)[1]
  else if(a > max(AA)) xi.a <- Xi[B,] #jitter(Xi)[B]
  else {
      
    r <- sum(AA < a)
    w <- (a - AA[r]) / (AA[r+1] - AA[r])
    xi.a <- (1 - w) * Xi[r+1,] + w * Xi[r,]
      
  }
  V <- eJ$vectors %*% diag(sqrt(xi.a / eJ$values))
  return(mle + sqrt(qchisq(1 - a, 2)) * V %*% u)
    
}
start.time <- Sys.time()
numCores <- detectCores() 
clus <- makeCluster(numCores)
registerDoParallel(clus)
lTh <- t( foreach(k=1:N, .combine='cbind') %dopar% { H(k) } )
stopCluster(clus)
end.time <- Sys.time()
print(as.numeric(end.time - start.time))
  

# Plot to visualize the samples

plot(lTh[,1], lTh[,2], xlab=expression(paste("log ", theta[1])), 
       ylab=expression(paste("log ", theta[2])))
  

# Transform samples to a possibility contour

f <- function(z) logrellik(z, X)
f.lTh <- apply(lTh, 1, f)
lth1 <- seq(min(lTh[,1]), max(lTh[,1]), length=100)
lth2 <- seq(min(lTh[,2]), max(lTh[,2]), length=100)
approx.pl.lik <- outer(lth1, lth2)
for(i in seq_along(lth1)) {
    
  for(j in seq_along(lth2)) {
      
    v <- c(lth1[i], lth2[j])
    approx.pl.lik[i,j] <- mean(f.lTh <= f(v))
      
  }
    
}

plot(0, 0, type="n", xlim=c(4, 17), ylim=c(5, 31), 
    xlab=expression(theta[1]), ylab=expression(theta[2]))
contour(exp(lth1), exp(lth2), approx.pl.lik, drawlabels=FALSE, add=TRUE)

  
  

