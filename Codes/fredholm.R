
# Some R codes related to the paper 
#
# M. Chae, R. Martin, and S. G. Walker (2019). On an algorithm for solving Fredholm equations 
# of the first kind. Statistics and Computing, volume 29, pages 645-654; arXiv:1709.02695


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

newpr <- function(x, fx, u, K, p0, M) {

  if(length(u) %% 2 != 0) stop("u-grid length should be of even length!")
  if(length(x) %% 2 != 0) stop("x-grid length should be of even length!")
  tK <- t(K)
  p <- p0
  Kp <- t(tK * p)
  fp <- apply(Kp, 1, myint, x=u)
  D <- myint(fx * log(fx / fp), x=x)
  #D.diff <- numeric(M)
  for(m in 1:M) {

    p <- apply(Kp * fx / fp, 2, myint, x=x)
    Kp <- t(tK * p)
    fp <- apply(Kp, 1, myint, x=u)
    Dm <- myint(fx * log(fx / fp), x=x)
    #D.diff[m] <- D[length(D)] - Dm
    D <- c(D, Dm)

  }
  return(list(u=u, pu=p, x=x, fpx=fp, D=D))

}


## Galaxy data illustration

M <- 25
library(MASS)
data(galaxies)
Y <- galaxies / 1000
d <- density(Y, n=2**9)
y <- x <- d$x
fy <- d$y
p0 <- dt(x, df=1)
K <- dnorm(outer(y, x, "-"))
o <- newpr(y, fy, x, K, p0, M)
plot(x, o$pu, type="l", xlab=expression(theta), ylab=expression(p(theta)))
ylim <- c(0, max(c(o$fpx, fy)))
hist(Y, freq=FALSE, breaks=10, xlab="x", ylab="f(x)", col="gray", border="white", main="", ylim=ylim)
lines(y, o$fpx, lty=2)
lines(y, fy)
plot(o$D, type="l", xlab="m", ylab=expression(D(f,f[m])))


