
#--------------------------------------------------------------------------#
# R code for computing optimal plausibility function for Poisson mean.     #
# Authors: Ryan Martin, Duncan Ermini Leaf, and Chuanhai Liu.              #
# For details, see our paper                                               #
#                                                                          #
#  "Optimal inferential models for a Poisson mean," 2012.                  #
#                                                                          #
# Submitted for publication.                                               #
#--------------------------------------------------------------------------#


pl.pois <- function (theta, x, PRS=c("Default0", "Default1", "Score-Balanced"), b=0) {

  sb.order <- function(X, lambda, d) {

    x <- X
    if (lambda == 0) return(structure(x + 1, names = x))
    if(missing(d)) d <- dpois(x, lambda)
    r <- integer(length(x))
    dd <- x - lambda
    ES <- SD <- SV <- 0
    for (i in 1:length(r)) {

      if (all(dd > 0) || all(dd < 0)) {

        r[i:length(r)] <- X[order(abs(dd))] + 1
        break

      }
      ddp <- dd
      ddp[dd < 0] <- max(dd)
      ddn <- dd
      ddn[dd > 0] <- min(dd)
      K <- c(which.min(ddp), which.max(ddn))
      V <- SV + ((X[K] - lambda)**2 - lambda) * d[K]
      if(any(V > 0)) {

        if(all(V > 0)) stop("Whoops: something's wrong!") else K[V > 0] <- K[V < 0]

      }
      z <- (ES + d[K] * dd[K]) / (SD + d[K])
      k <- K[which.min(abs(z))]
      ES <- ES + d[k] * dd[k]
      SD <- SD + d[k]
      SV <- SV + ((X[k] - lambda)**2 - lambda) * d[k]
      r[i] <- X[k] + 1
      X <- X[-k]
      dd <- dd[-k]
      d <- d[-k]

    }
    if (length(unique(r)) != length(x)) stop("Whoops: something's wrong!")
    names(r) <- x[r]
    return(r)

  }

  pl.default <- function (theta, x) {

    pl <- numeric(length(theta))
    for(i in seq_along(theta)) {

      pl.1 <- max(1 - 2 * pgamma(theta[i], x), 0)
      pl.2 <- max(2 * pgamma(theta[i], x + 1) - 1, 0)
      pl[i] <- 1 - pl.1 - pl.2

    }
    pl <- data.frame(theta=theta, pl=pl)
    attr(pl, "x") <- x
    return(pl)

  }

  default.order <- function(X, lambda, d) {

    p <- ppois(X, lambda = lambda)
    k <- which(p >= 0.5)[1]
    if(k >= length(p)) stop("Whoops: k >= length(X)")
    s <- (k+1):length(p)
    p[s] <- 1 - p[s-1]
    r <- order(p, decreasing=TRUE)
    names(r) <- X[r]
    return(r)

  }

  im.order <- if(PRS=="Score-Balanced") sb.order else
              if(PRS=="Default0") return(pl.default(theta, x)) else
              if(PRS=="Default1") default.order else
              stop("unknown PRS", PRS)

  pl <- numeric(length(theta))
  for(i in seq_along(theta)) {

    lambda <- theta[i]
    d <- dpois(X <- 0:(2*max(as.integer(qpois(0.999999, lambda)), x+1)), lambda)
    r <- im.order(X, lambda, d)
    names(d) <- X
    a <- d[r]
    a <- c("-1"=0,a)
    k <- which(names(a)==x)
    pl[i] <- 1 - (bel <- sum(a[1:(k-1)])) # cumsum(a)[k-1]
    # EB-SB effect: put all conflict mass on lambda=b
    # if lambda is on the boundary and there are any conflict cases
    if(lambda == b && ppois(X[r[1]] - 1, lambda) >= ppois(x, b)) pl[i] = 1.0

  }
  pl <- data.frame(theta=theta, pl=pl)
  attr(pl, "x") <- x
  if(length(theta) == 1) attr(pl, "r") <- as.integer(names(r))
  return(pl)

}



