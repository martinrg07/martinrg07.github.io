
# Authors: Ryan Martin (martinrg@purdue.edu) and Sahil Patel (spate27@ncsu.edu)
# Reference: Probabilistic approximations of possibility measures 
#            and fiducial inference for discrete models (arXiv:...)

# This code is focused on extracting the inner probabilistic approximation of
# the possibilistic inferential model (IM), as defined in the paper, specifically 
# for the binomial model, which is discrete.  This involves first finding the  
# jump points in the IM's contour function, as in the get.jumps() function below, 
# which are the candidate atoms in the inner probabilistic approximation.  The 
# rest of the steps can be found in the get.ipa.binom() function. 

# Auxiliary functions 

expit <- function(u) exp(u) / (1 + exp(u))

int <- function(f, x, tol=1e-10) {
  
  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) out <- sum((f[-1] + f[-n]) * (x[-1] - x[-n])) / 2
  else out <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n-1)]) + 2 * sum(f[seq(2, n-1, 2)]) )
  return(out)
  
}


# Extract the jumps in the binomial contour

get.jumps.binom <- function(n, x) {
  
  mle.x <- x / n
  g <- function(u) {
    
    o <- 0 * u
    uu <- u[u > 0 & u < 1]
    o[u > 0 & u < 1] <- uu * log(uu) + (1 - uu) * log(1 - uu)
    return(o)
    
  }
  if(x == 0) y <- 1:n
  else if(x == n) y <- 0:(n - 1)
  else y <- c(0:(x - 1), (x + 1):n)
  mle.y <- y / n
  G <- (g(mle.y) - g(mle.x)) / (mle.y - mle.x)
  jmp <- expit(G)
  return(sort(jmp))
  
}


# Binomial log-relative likelihood

binom.lrl <- function(y, n, p0) {
  
  q0 <- 1 - p0
  if(y == 0) t1 <- 0 else t1 <- y * (log(n * p0) - log(y))
  if(y == n) t2 <- 0 else t2 <- (n - y) * (log(n * q0) - log(n - y))
  return(t1 + t2)
  
}


# Binomial possibility contour 

binom.pl <- function(y, n, th) {
  
  s <- 0:n
  f <- function(x) binom.lrl(x, n, th)
  log.t <- f(y)
  log.T <- sapply(s, f)
  return(sum(dbinom(s[log.T <= log.t], n, th)))
  
}

# Exract the inner probabilistic approximation
# Note: The present code does not directly handle the case of a purely symmetric 
#       contour, i.e., where x = n / 2.  Code that handles the symmetric case 
#       will be added later. 

get.ipa.binom <- function(n, x, th) {
  
  pl <- function(z) binom.pl(x, n, z)
  
  eps <- 0.001
  mle <- x / n
  jmps <- get.jumps.binom(n, x)
  
  pl.lo <- sapply(jmps + eps * sign(jmps - mle), pl)
  pl.hi <- sapply(jmps - eps * sign(jmps - mle), pl)
  
  ix <- order(pl.hi, pl.lo)
  J <- length(pl.hi)
  
  q.disc <- 0 * jmps
  
  for(j in 1:(J-2)) {
    
    i <- ix[j]
    o <- (pl.lo < pl.hi[i]) & (pl.hi >= pl.hi[i])
    if(any(o)) q.disc[i] <- pl.hi[i] - max(pl.lo[o])
    
  }
  q.disc[ix[J]] <- q.disc[ix[J-1]] <- (1 - max(pl.lo)) / 2
  
  dmass <- sum(q.disc)
  
  f <- function(z) as.numeric(any((z >= pl.lo) & (z <= pl.hi)))
  
  dpl <- function(th) {
    
    s <- 0:n
    g <- function(y) binom.lrl(y, n, th)
    log.t <- g(x)
    log.T <- sapply(s, g)
    ss <- s[log.T <= log.t]
    score <- ss / th - (n - ss) / (1 - th)
    return(sum(score * dbinom(ss, n, th)))
    
  }
  
  pl.th0 <- sapply(th, pl)
  pl.th <- pl.th0
  pl.th[th <= mle] <- cummax(pl.th[th <= mle])
  pl.th[th > mle] <- rev(cummax(rev(pl.th[th > mle])))
  
  dpl.th <- sapply(th, dpl)
  
  tmp <- ifelse(th <= x / n, pmax(dpl.th, 0), -pmin(dpl.th, 0))
  q.cont <- 0.5 * (1 + sapply(pl.th, f)) * tmp
  norm <- int(q.cont, th)
  
  out.cont <- cbind(th=th, pdf=q.cont / norm)
  out.disc <- cbind(atoms=jmps, wt=q.disc, pmf=q.disc / dmass)
  out <- list(th=th, pl=pl, cont=out.cont, disc=out.disc)
  return(out)
  
}


# Illustrative example

n <- 40
x <- 0.6 * n
  
th <- seq(0.01, 0.99, length=257)
pl.th <- sapply(th, function(z) binom.pl(x, n, z))
plot(th, pl.th, type="l", xlab=expression(theta), ylab=expression(pi[x](theta)))
  
o <- get.ipa.binom(n, x, th)
plot(o$disc[,1], o$disc[,3], type="h", xlab=expression(theta), ylab="Mass")
plot(o$cont[,1], o$cont[,2], type="l", xlab=expression(theta), ylab="Density")
  
d.mass <- sum(o$disc[,2])
  
cdf.disc <- sapply(th[-c(1,2)], function(z) sum(o$disc[o$disc[,1] <= z, 3]))
cdf.cont <- sapply(th[-c(1,2)], function(z) int(o$cont[o$cont[,1] <= z, 2], o$cont[o$cont[,1] <= z, 1]))
cdf <- d.mass* cdf.disc + (1 - d.mass) * cdf.cont
plot(th[-c(1,2)], cdf, type="l", xlab="th", ylab="CDF")
mle <- x / n
sd.mle <- sqrt(mle * (1 - mle) / n)
curve(pnorm(x, mle, sd.mle), add=TRUE, lty=2)


# Simulation study
  
reps <- 1000
n <- 40
theta <- 0.6
X <- rbinom(reps, n, theta)
th <- seq(0.01, 0.99, length=157)
out <- numeric(reps)
for(r in 1:reps) {
    
  x <- X[r]
  R.theta <- binom.lrl(x, n, theta)
  o <- get.ipa.binom(n, x, th)
  atoms <- o$disc[,1]
  mass <- o$disc[,3]
  R.atoms <- binom.lrl(x, n, atoms)
  out[r] <- sum((R.atoms <= R.theta) * mass)
    
}
plot(sort(out), (1:reps) / reps, type="l", main="", xlab="u", ylab="CDF of U")
abline(a=0, b=1, lty=3)





