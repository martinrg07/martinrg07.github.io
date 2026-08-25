

## R code to implement the recursive Dirichlet process mixture of normals approximation 
## to the predictive distribution, proposed in 
##
##   Hahn, Martin, and Walker (201x). On recursive Bayesian predictive distributions, 
##   Journal of the American Statistical Association, to appear; arXiv:1508.07448
##
## The specific insuracne examples here are presented in the paper 
##
##   Hong and Martin (201x). Real-time Bayesian nonparametric prediction of solvency 
##   risk, submitted for publication; SSRN Abstract 2984878
##
## Examples below use two standard insurance datasets (Danish and Norwegian fire insurance
## claims) that are available in the accompanying "dpmrec.Rd" file.


## Auxiliary function for numerical integration

myint <- function(fx, x, m=2) {
  
  h <- x[2] - x[1]
  n <- length(x)
  if(!(m %in% c(0, 1, 2))) stop("Only Riemann, Trapezoid, or Simpson is implemented!")
  if(m == 0) intf <- h * sum(fx[-n])
  else if(m == 1) intf <- (h / 2) * sum(fx[-1] + fx[-n])
  else {
    
    if(n %% 2 != 0) stop("Simpson's rule requires grid of even length!")
    intf <- (h / 3) * (sum(fx) + sum(fx[2:(n-1)]) + 2 * sum(fx[seq(2, n-1, by=2)]))
    
  }
  return(intf)
  
}


## Main function: recursive DP mixture (of normals)
## INPUT:
## X = data sequence
## xgrid = points on which the predictive distribution is to be estimated
## P0 = starting guess
## rho = correlation parameter
## q = quantile, if VaR and CTE calculations are desired
## do.cte = logical, if FALSE, don't do CTE calculations

## OUTPUT:
## xgrid = same as xgrid input
## P = final estimated distribution function
## p = final estimated density function
## var = value-at-risk trajectory, depending on q
## cte = conditional tail expectation trajectory, depending on q

dpmrec <- function(X, xgrid, P0, rho=0.9, q=NULL, do.cte=FALSE) {
  
  n <- length(X)
  alpha <- 1 / (1 + 1:n)
  s <- sqrt(1 - rho**2)
  if(is.function(P0)) F <- P0(xgrid) else P <- P0
  if(is.null(q)) vr <- cte <- ix <- NULL else {
    
    vr <- rep(NA, n)
    lo <- mean(xgrid)
    hi <- max(xgrid)
    if(do.cte) { cte <- rep(NA, n); xx <- xgrid[-1] }
    
  }
  for(i in 1:n) {
    
    P.fun <- approxfun(xgrid, P)
    if(!is.null(q)) { #
      
      g <- function(x) P.fun(x) / max(P) - q
      vr[i] <- uniroot(g, lower=lo, upper=hi, extendInt="yes")$root
      if(do.cte && (i > round(0.1 * n))) {
        
        p <- diff(P) / diff(xgrid)
        ix <- xx > vr[i]
        if(sum(ix) %% 2 == 0) m <- 2 else m <- 1
        cte[i] <- myint((xx * p)[ix], xx[ix], m) / (1 - q)
        
      }
      
    }
    Pn.X <- P.fun(X[i])
    P <- (1 - alpha[i]) * P + alpha[i] * pnorm((qnorm(P) - rho * qnorm(Pn.X)) / s)
    
  }
  return(list(xgrid=xgrid, P=P, p=diff(P) / diff(xgrid), var=vr, cte=cte))
  
}


if(FALSE) {
  
  load("dpmrec.Rd")
  
  # Example 1: Norwegian 1988 fire claims data
  
  norway <- dpmrec.data$norway
  loss <- norway[, 1]
  year <- norway[, 2]
  X <- log(loss[year == 88]) 
  set.seed(77)
  X <- sample(X)
  xgrid <- seq(-15, 15, len=1000) 
  norway.out.1 <- dpmrec(X, xgrid, pt(xgrid, df=2), 0.9, 0.995, TRUE)
  norway.out.2 <- dpmrec(X, xgrid, pnorm(xgrid, 0, 4), 0.9, 0.995, TRUE)

  ylim <- range(c(norway.out.1$var, norway.out.2$var))
  plot(norway.out.1$var, xlab="Index", ylab="VaR", type="l", ylim=ylim)
  lines(norway.out.2$var, lty=2)
  
  ylim <- c(0, max(norway.out.1$p, norway.out.2$p))
  hist(X, freq=FALSE, breaks=30, col="gray", border="white", main="", xlab="log(loss)", ylim=ylim)
  lines(norway.out.1$xgrid[-1], norway.out.1$p)
  lines(norway.out.2$xgrid[-1], norway.out.2$p, lty=2)
  
  
  # Example 2: Danish fire insurance data
  
  danish <- dpmrec.data$danish
  set.seed(77)
  X <- sample(log(danish)) 
  xgrid <- seq(-12, 12, len=1000) 
  danish.out.1 <- dpmrec(X, xgrid, pt(xgrid, df=2), 0.95, 0.995, TRUE)
  danish.out.2 <- dpmrec(X, xgrid, pnorm(xgrid, 0, 4), 0.95, 0.995, TRUE)
  
  ylim <- range(c(danish.out.1$var, danish.out.2$var))
  plot(danish.out.1$var, ylab="VaR", type="l", ylim=ylim)
  lines(danish.out.2$var, lty=2)
  
  ylim <- c(0, max(danish.out.1$p, danish.out.2$p))
  hist(X[X <= log(50)], freq=FALSE, breaks=45, col="gray", border="white", main="", xlab="log(loss)", ylim=ylim)
  lines(danish.out.1$xgrid[-1], danish.out.1$p)
  lines(danish.out.2$xgrid[-1], danish.out.2$p, lty=2)
  
}







