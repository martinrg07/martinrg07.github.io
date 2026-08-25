
##--------------------------------------------------------##
## R source code for examples in the "IM Basics" paper    ##
## Authors: Ryan Martin and Chuanhai Liu                  ##
##--------------------------------------------------------##


# Bisection method

bisection <- function (f, a, b, eps=1e-08, maxiter=1000, verbose=FALSE, ...) {

  x <- (a + b) / 2
  t <- 0
  if(verbose) cat("\n")

  repeat {

    t <- t + 1
    if(f(a, ...) * f(x, ...) <= 0) b <- x else a <- x
    x.new <- (a + b) / 2
    if(verbose) cat("[", t, "]", x.new, f(x.new, ...), "\n")
    if(abs(x.new - x) < eps | t >= maxiter) {

      if(t >= maxiter) warning("Maximum number of iterations reached!")
      break

    }
    x <- x.new

  }

  if(verbose) cat("\n")
  out <- list(solution=x.new, value=f(x.new, ...), iter=t)
  return(out)

}


# Plot of the plausibility function in Poisson example

pois.plot <- function(x=5, tt=NULL) {

  T <- seq(0, 13, len=300)
  P <- 0 * T
  for(j in 1:length(T)) {

    G1 <- pgamma(T[j], shape=x+1, scale=1)
    G2 <- pgamma(T[j], shape=x+0, scale=1)
    P[j] <- 1 - max(1 - 2 * G2, 0) - max(2 * G1 - 1, 0)

  }
  plot(x=0, y=0, type="n", xlim=range(T), ylim=c(0,1), xlab=expression(theta), ylab=expression(pl[x](theta)))
  lines(T, P, lwd=2)
  if(!is.null(tt)) {

    pl.tt <- 0 * tt
    pval.tt <- 0 * tt
    for(j in 1:length(tt)) {

      G1 <- pgamma(tt[j], shape=x+1, scale=1)
      G2 <- pgamma(tt[j], shape=x+0, scale=1)
      pl.tt[j] <- 1 - max(1 - 2 * G2, 0) - max(2 * G1 - 1, 0)
      pval.tt[j] <- 1 - dpois(x, lambda=tt[j])

    }
    print(cbind(t=tt, pl.t=pl.tt, pval.t=pval.tt))

  }

}


# Plausibility plot for the Gaussian example

gauss.plot <- function(x=5, tt=NULL) {

  T <- seq(x - 3, x + 3, len=300)
  pl <- 1 - abs(2 * pnorm(x - T) - 1)
  plot(T, pl, ylim=c(0, 1), type="l", lwd=2, xlab=expression(theta), ylab=expression(pl[x](theta)))

}


# Plausibility interval for the Poisson example

pois.pi <- function(x=5, alpha=0.10) {

  pl <- function(t) 1 - max(1 - 2 * pgamma(t, x), 0) - max(2 * pgamma(t, x+1) - 1, 0)
  f <- function(t) pl(t) - alpha
  pi.lower <- bisection(f=f, a=0, b=x)$solution
  pi.upper <- bisection(f=f, a=x, b=3*x)$solution
  return(c(pi.lower, pi.upper))

}


# "Optimal" plausibility function in exponential example

exp.score.bal.int <- function(x, theta) {

  f <- function(xx) xx * exp(-xx / theta) - x * exp(-x / theta)
  if(x < theta) {

    a <- x
    b <- bisection(f, theta, qexp(0.99, 1 / theta))$solution

  } else {

    b <- x
    a <- bisection(f, 0, theta)$solution

  }
  return(c(a, b))

}


exp.demo <- function(tmax=2, theta0=1) {

  tt <- seq(-1, tmax, len=50)
  V <- function(t) t**2 - t / 2 / theta0 - 1 / theta0**2
  h <- function(t) theta0 * dexp(theta0 * t + theta0)
  th <- function(t) t * h(t)
  B.upper <- B.lower <- seq(0, max(V(tt)), by=max(V(tt)) / 50)
  plot(tt, V(tt), col="white", xlab="t", ylab="V(t)")
  for(i in seq_along(B.lower)) {

    B.lower[i] <- (exp.score.bal.int(theta0 * B.upper[i] + theta0, theta0)[1] - theta0) / theta0
    ix <- integrate(th, 0, 0.5)$value + integrate(th, B.lower[i], 0)$value > 0
    lines(x=c(B.lower[i], B.upper[i]), y=rep(B.upper[i], 2), col=ifelse(ix, "black", "gray"))

  }
  lines(tt, V(tt) * h(tt), col="gray")
  lines(tt, tt * h(tt))
  lines(tt, V(tt), lwd=2)
  abline(h=0)
  abline(v=c(0, 0.5), h=-1, lty=3)

}


exp.pl.opt <- function(x, theta) {

  pl <- 0 * theta
  for(i in seq_along(theta)) {

    bal.int <- exp.score.bal.int(x, theta[i])
    pl[i] <- 1 - pexp(bal.int[2], 1 / theta[i]) + pexp(bal.int[1], 1 / theta[i])

  }
  pl.def <- 1 - abs(2 * pexp(x, 1 / theta) - 1)
  plot(theta, pl.def, type="l", lwd=2, col="gray", xlab=expression(theta), ylab=expression(pl[x](theta)))
  lines(theta, pl, lwd=2)
  abline(h=0.1, v=x, lty=3, lwd=2)
  return(data.frame(theta=theta, pl.def=pl.def, pl=pl))

}



# "Optimal plausibility function for the Poisson example

pois.pl.opt <- function (theta, x, PRS=c("Default0", "Default1", "Score-Balanced")) {

    sb.order <- function(X, lambda, d) {
        x <- X
        if (lambda == 0)
            return(structure(x + 1, names = x))
        if(missing(d)) d <- dpois(x, lambda)
        r <- integer(length(x))
        dd <- x/lambda - 1
        ES <- 0
        SD <- 0
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
            z <- (ES + d[K] * dd[K])/(SD + d[K])
            k <- K[which.min(abs(z))]
            ES <- ES + d[k] * dd[k]
            SD <- SD + d[k]
            r[i] <- X[k] + 1
            X <- X[-k]
            dd <- dd[-k]
            d <- d[-k]
        }
        if (length(unique(r)) != length(x))
            stop("Whoops: somethings is wrong")
        names(r) <- x[r]
        r
    }

  pl.default <- function (theta, x) {
    pl <- numeric(length(theta))
    for(i in seq_along(theta))
        pl[i] <- (1 - max(1 - 2 * pgamma(theta[i], x), 0) - max(2 * pgamma(theta[i], x + 1) - 1, 0))
    pl <- data.frame(theta=theta, pl = pl)
    attr(pl, "x") <- x
    pl
  }

  default.order <- function(X, lambda, d){
     p <- ppois(X, lambda = lambda)
     k <- which(p>=0.5)[1]
     if(k>=length(p)) stop("Whoops: k >= length(X)")
     s <- (k+1):length(p)
     p[s] <- 1 - p[s-1]
     r <- order(p, decreasing=TRUE)
     names(r) <- X[r]
     r
  }

  im.order <- if((PRS <- PRS[1])=="Score-Balanced")  sb.order
                else if (PRS=="Default0") return(pl.default(theta, x))
                else if (PRS=="Default1") default.order
                else stop("unknown PRS", PRS)

  pl <- numeric(length(theta))
  for(i in seq_along(theta)){

    lambda <- theta[i]
    d <- dpois(X <- 0: max(as.integer(qpois(0.999999, lambda)), x+1), lambda)
    r <- im.order(X, lambda, d)
    names(d) <- X
    a <- d[r]
    a <- c("-1"=0,a)
    k <- which(names(a)==x)
    pl[i] <- 1 - (bel <- sum(a[1:(k-1)])) # cumsum(a)[k-1]

  }
  pl <- data.frame(theta=theta, pl = pl)
  attr(pl, "x") <- x
  return(pl)

}

pois.pl.opt.plot <- function(x, theta) {

  pl.def <- pois.pl.opt(theta, x, "Default0")
  pl.sb <- pois.pl.opt(theta, x, "Score-Balanced")
  plot(pl.def, type="l", lwd=2, col="gray", xlab=expression(theta), ylab=expression(pl[x](theta)))
  lines(pl.sb, lwd=2)
  abline(h=0.1, lwd=2, lty=3)

}



# Simulation for the "many-exponential-rates" problem.

mer.sim <- function(n=c(10,90), thetas=1:5, alpha=0.05, num.simu=5000, mcmc.size=50000) {

   N <- sum(n)
   M <- N - 1
   b <- 1 / ((1:M) - 1 + 0.7); a <- rev(b)
   aa <- bb <- 0.1
   UU <- matrix(rexp(mcmc.size * N), mcmc.size, N)

   f <- function(u) {

     v <- cumsum(u)
     v <- v[-length(v)] / v[length(v)]
     return( sum(a * log(v) + b * log(1-v)) )

   }

   lr <- function(u) sum(log(u)) - length(u) * log(mean(u))

   bf <- function(u) {

      m0 <- exp(aa * log(bb) + lgamma(aa + sum(n)) - (aa + sum(n)) * log(sum(u) + bb) - lgamma(aa))
      m1.0 <- function(t) {

        lp1 <- sum(n) * log(t) + dgamma(t, shape=aa, rate=bb, log=TRUE)
        lp2 <- apply(log(outer(t, u, "+")), 1, sum)
        return(exp(lp1 - 2 * lp2))

      }
      m1 <- integrate(m1.0, lower=0, upper=Inf)$value
      return(m1 / m0)

   }

   F <- sort(apply(UU, 1, f))
   L <- sort(apply(UU, 1, lr))
   power <- data.frame(theta=thetas, im=0 * thetas, lr=0 * thetas)
   k <- 0

   for(theta in thetas){

     p.im <- p.lr <- numeric(num.simu)

     for(i in 1:num.simu) {

       x <- c(rexp(n[1]), rexp(n[2]) / theta)
       fx <- f(x)
       lrx <- lr(x)
       p.im[i] <- mean(fx >= F)
       p.lr[i] <- mean(lrx >= L)
       cat("\r", theta, i)

     }

     k <- k + 1
     power$im[k] <- mean(p.im <= alpha)
     power$lr[k] <- mean(p.lr <= alpha)

   }

   plot(x=thetas, y=power$im, xlab=expression(theta), ylab="Power", ylim=c(0,1), type="b", lwd=2, pch=1)
   lines(x=thetas, y=power$lr, type="b", lwd=2, pch=4)
   if(draw.leg) legend(x="bottomright", inset=0.05, lwd=2, pch=c(1,2,4), c("IM new", "IM old", "LR"))
   cat("\n")
   return(list(n=n, alpha=alpha, power=power))

}
