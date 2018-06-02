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

  # To convert mu, sigma, nu and tau to original parameters
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

#---------------------------------------------------------------------

integrate(dWGEE, lower=0, upper=Inf,
          mu=5, sigma=0.5, nu=1)

# Figure 3
curve(dWGEE(x, mu=5, sigma=0.5, nu=1),
      from=0, to=6, ylab='f(x)', las=1)
curve(dWGEE(x, mu=1, sigma=0.5, nu=1),
      from=0, to=6, ylab='f(x)', las=1, add=TRUE)
curve(dWGEE(x, mu=0.1, sigma=0.5, nu=1),
      from=0, to=6, ylab='f(x)', las=1, add=TRUE)

# The incomplete beta function
ibeta <- function(x, a, b) {pbeta(x, a, b, lower.tail=F) * beta(a, b)}

pWGEE <- function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  # To convert mu, sigma, nu and tau to original parameters
  alp <- mu
  bet <- sigma
  lam <- nu
  
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

curve(pWGEE(x, mu=5, sigma=0.5, nu=1),
      from=0, to=6, ylab='F(x)', las=1)

# Test
integrate(dWGEE, lower=0, upper=1,
          mu=5, sigma=0.5, nu=1)

pWGEE(q=1, mu=5, sigma=0.5, nu=1)



