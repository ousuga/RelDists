# Weighted Generalized Exponential-Exponential distribution
# Taken from http://www.ijmex.com/index.php/ijmex/article/view/267/188

dWGEE <- function(x, mu, sigma, nu, log=FALSE) {
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))

  # To convert mu, sigma, nu to original parameters
  alp <- mu
  bet <- sigma
  lam <- nu

  loglik <- log(bet * lam) - lam * x + (bet - 1) * log(1-exp(-lam*x)) +
    log(1 - exp(-lam * alp * x)) - log(1 - bet * beta(alp + 1, bet))
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}

# Examples

# The area under the pdf
integrate(dWGEE, lower=0, upper=Inf, mu=5, sigma=0.5, nu=1)
integrate(dWGEE, lower=0, upper=Inf, mu=1, sigma=0.5, nu=1)
integrate(dWGEE, lower=0, upper=Inf, mu=0.1, sigma=0.5, nu=1)

# Figure 3.a
curve(dWGEE(x, mu=5, sigma=0.5, nu=1),
      from=0, to=6, ylab='f(x)', las=1)
curve(dWGEE(x, mu=1, sigma=0.5, nu=1), col='tomato', add=TRUE)
curve(dWGEE(x, mu=0.1, sigma=0.5, nu=1), col='blue', add=TRUE)

#---------------------------------------------------------------------
pWGEE <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  # To convert mu, sigma, nu to original parameters
  alp <- mu
  bet <- sigma
  lam <- nu
  
  # The incomplete beta function
  ibeta <- function(x, a, b) {pbeta(x, a, b, lower.tail=F) * beta(a, b)}
  
  cdf <- ((1-exp(-lam*q))^bet - bet * ibeta(exp(-lam*q), alp+1, bet)) / 
    (1 - bet * beta(alp + 1, bet))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf 
}

# Plotting 3 F(x)
curve(pWGEE(x, mu=5, sigma=0.5, nu=1),
      from=0, to=6, ylab='F(x)', las=1)
curve(pWGEE(x, mu=1, sigma=0.5, nu=1), col='tomato', add=TRUE)
curve(pWGEE(x, mu=0.1, sigma=0.5, nu=1), col='blue', add=TRUE)


#---------------------------------------------------------------------
qWGEE <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  F.inv <- function(y, mu, sigma, nu) {
    uniroot(function(x) {pWGEE(x,mu,sigma,nu) - y},
            interval=c(0, 99999))$root
  }
  F.inv <- Vectorize(F.inv)
  F.inv(p, mu, sigma, nu)
}


rWGEE <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  # To convert mu, sigma, nu to original parameters
  alp <- mu
  bet <- sigma
  lam <- nu
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qWGEE(p, mu, sigma, nu)
  r
}

rWGEE(n=10, mu=5, sigma=0.5, nu=1)

hist(rWGEE(n=1000, mu=5, sigma=0.5, nu=1), las=1, freq=FALSE)
curve(dWGEE(x, mu=5, sigma=0.5, nu=1), from=0, to=8,
      col='purple', add=TRUE, lwd=2)

