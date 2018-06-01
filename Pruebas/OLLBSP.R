
dOLLBSP <- function(x, mu, sigma, nu, tau, log=FALSE) {
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must be positive", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must be positive", "\n", ""))
  # To convert mu, sigma, nu and tau to original parameters
  alp <- mu
  bet <- sigma
  a   <- nu
  b   <- tau
  # An auxiliar variable v 
  v <-  a^-1 * (sqrt(x/b) - 1/sqrt(x/b))
  
  loglik <- log(alp) + log(bet) - (log(2*a) + a^2 + 0.5 * log(2*pi*b)) - 
    1.5 * log(x) + log(x+b) - (x/b + b/x)/(2*a^2) +
    (alp-1) * log(pnorm(v)) + (alp-1) * log(1-pnorm(v)) +
    bet * pnorm(v)^alp / (pnorm(v)^alp + (1-pnorm(v))^alp) -
    log(exp(bet) - 1) - 2 * log(pnorm(v)^alp + (1-pnorm(v))^alp)
  if (log == FALSE) 
    density <- exp(loglik)
  else 
    density <- loglik
  return(density)
}

curve(dOLLBSP(x, mu=0.25, sigma=0.5, nu=0.1, tau=0.5),
      from=0.3, to=1, ylab='f(x)', las=1)

integrate(dOLLBSP, lower=0.3, upper=0.8,
          mu=0.25, sigma=0.5, nu=0.1, tau=0.5)

curve(dOLLBSP(x, mu=1.5, sigma=0.1, nu=0.1, tau=0.5),
      from=0.4, to=0.6, ylab='f(x)', las=1) 

curve(dOLLBSP(x, mu=0.05, sigma=3.5, nu=0.1, tau=0.5),
      from=0.1, to=1, ylab='f(x)', las=1) 


