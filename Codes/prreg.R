
##-----------------------------------------------------------------------------------
## TITLE:      R codes for hybrid predictive recursion-EM algorithm for regression
##
## AUTHORS:    R. Martin (rgmartin@uic.edu) and Z. Han
##
## REFERENCE:  Martin and Han, "Robust regression via predictive recursion maximum
##             likelihood", arXiv:1306.3185
##
## NOTE:       A version that implements the PR portion of PR-EM in C would run
##             faster and may be created later.
##-----------------------------------------------------------------------------------


library(MASS)


int <- function(f, x, tol=1e-10) {

  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) o <- (1 / 2) * sum((f[-1] + f[-n]) * (x[-1] - x[-n]))
  else o <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n - 1)]) + 2 * sum(f[seq(2, n - 1, 2)]) )
  return(o)

}


pr.reg <- function(X, d, U, f0, w, perm) {

  n <- length(X)
  if(missing(f0)) { f0 <- 1 + 0 * U; f0 <- f0 / int(f0, U) }
  if(missing(w)) w <- function(i) 1 / (i + 1)
  N <- ncol(perm)
  f.avg <- 0 * U
  L.avg <- 0
  wt <- 0 * X
  for(j in 1:N) {

    f <- f0
    L <- 0
    for(i in 1:n) {

      num <- d(X[perm[i, j]], U) * f
      den <- int(num, U)
      wt[perm[i, j]] <- wt[perm[i, j]] + int(num / U**2, U) / den
      L <- L + log(den)
      f <- (1 - w(i)) * f + w(i) * num / den

    }
    f.avg <- f.avg + f / N
    L.avg <- L.avg + L / N

  }
  return(list(f=f.avg, L=-L.avg, wt=wt / N))

}



prem.reg <- function(X, Y, perm, get.ci=FALSE) {

  tol <- 1e-05
  ols <- lm(Y ~ X - 1)
  sigma <- summary(ols)$sigma
  UUmax <- max(5 * sigma, 50)
  UU <- seq(1e-05, UUmax, len=201)
  K <- function(z, u) dnorm(z, 0, u)
  if(missing(perm)) {

    n <- length(Y)
    nperm <- 25
    perm <- matrix(0, ncol=nperm, nrow=n)
    perm[,1] <- 1:n
    for(j in 2:nperm) perm[,j] <- sample(n)

  }
  else nperm <- ncol(perm)
  b.old <- as.numeric(ols$coefficients)
  f <- 1 + 0 * UU
  f <- f / int(f, UU)
  lik <- c()
  repeat {

    Z <- as.numeric(Y - X %*% b.old)
    o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
    lik <- c(lik, -o$L)
    ww <- o$wt
    b <- as.numeric(lm(Y ~ X - 1, weights=ww)$coefficients)
    if(sum(abs(b - b.old)) < tol) break else b.old <- b

  }
  Z <- as.numeric(Y - X %*% b)
  o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
  lik <- c(lik, -o$L)
  ww <- o$wt
  if(get.ci) {

    if(!exists("hessian")) library(numDeriv)   #requires the R package "numDeriv"
    prml <- function(bb) {

      Z <- as.numeric(Y - X %*% bb)
      L <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)$L
      return(L)

    }
    J <- solve(hessian(prml, b))
    ci <- b + 1.96 * outer(sqrt(diag(J)), c(-1, 1))

  } else ci <- NULL
  return(list(b=b, Yhat=as.numeric(X %*% b), ci=ci, f=o$f, UU=UU, wt=ww, lik=lik))

}


# Belgian phone call data in Example 1 of the paper

phone.example <- function() {

  library(MASS)
  data(phones)
  Y <- phones$calls
  X <- phones$year + 1900
  o <- prem.reg(X=cbind(1 + 0 * Y, X), Y=Y)
  plot(x=X, y=Y, xlab="Year", ylab="Calls (in millions)")
  abline(a=o$b[1], b=o$b[2])
  abline(lm(Y ~ X), lty=2)
  legend(x="topleft", inset=0.05, lty=1:2, c("PREM", "LS"))
  return(list(b=o$b, f=o$f, UU=o$UU, w=o$wt, lik=o$lik))

}


# TO RUN
# phone.example()











