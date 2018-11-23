#' @name GPW
#' 
#' @title 
#' The Generalized Power Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function , 
#' random generation and hazard function for the generalized power weibull distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu parameter one.    
#' @param sigma parameter two.
#' @param nu parameter three.        
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The generalized power weibull with parameters \code{mu}, \code{sigma}
#' and \code{nu} has density given by
#' 
#' f(x) = mu*sigma*nu^(-1)*x^(sigma-1)*(1+mu*(x^sigma))^(1/nu-1)*exp(1-(1+mu*(x^sigma))^(1/nu))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dGPW} gives the density, \code{pGPW} gives the distribution 
#' function, \code{qGPW} gives the quantile function, \code{rGPW}
#' generates random deviates and \code{hGPW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25), from = 0, to = 2.5, ylim = c(0, 3), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25), from = 0, to = 2.5, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25, lower.tail = FALSE), from = 0, to = 2.5, col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qGPW(p, mu = 0.5, sigma = 0.5, nu = 0.25), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rGPW(n = 10000, mu = 0.5, sigma = 0.5, nu = 0.25), freq = FALSE, xlab = "x", las = 1, main = "")
#' curve(dGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25),  from = 0, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hGPW(x, mu = 0.5, sigma = 0.5, nu = 0.25), from = 0, to = 6, ylim = c(0, 13), col = "red", las = 1, ylab = "The Hazard function")
GPW <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Generalized Power Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Power Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Generalized Power Weibull", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("GPW", "Generalized Power Weibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                 nopar = 3, 
                 type = "Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link    = as.character(substitute(nu.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun    = vstats$linkfun, 
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv    = vstats$linkinv, 
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr    = vstats$mu.eta, 
                 
                 
                 dldm = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dm <- y^sigma
                   dexp2dm <- (1/nu)*(exp1^((1/nu)-1))*dexp1dm
                   dldm <- 1/mu + ((1/nu)-1)*(dexp1dm/exp1) - dexp2dm
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dm <- y^sigma
                   d2exp2dm2 <- (1/nu)*((1/nu)-1)*(exp1^((1/nu)-2))*dexp1dm^2
                   d2ldm2 <- -(-(1/mu^2) - (((1/nu)-1)*dexp1dm^2)/exp1^2 - d2exp2dm2 )^2
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dd <- mu*(y^sigma)*log(y)
                   dexp2dd <- (1/nu)*(exp1^((1/nu)-1))*dexp1dd
                   dldd <- 1/sigma + log(y) + (((1/nu)-1)*(dexp1dd/exp1)) - dexp2dd
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dd <- mu*(y^sigma)*log(y)
                   d2exp1dd2 <- mu*(y^sigma)*(log(y))^2
                   d2exp2dd2 <- (1/nu)*( ((1/nu)-1)*(dexp1dd^2)*exp1^((1/nu)-2) + d2exp1dd2*exp1^((1/nu)-1))
                   d2ldd2 <- -(-(1/sigma^2) + (((1/nu)-1)/exp1^2)*(exp1*d2exp1dd2-dexp1dd^2) - d2exp2dd2)^2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp2dv <- -(1/nu^2)*(exp1^(1/nu))*log(exp1)
                   dldv <- -1/nu - log(exp1)/nu^2 - dexp2dv 
                   dldv 
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   d2exp2dv2 <- log(exp1)*( (2*exp1^(1/nu))/nu^3 + (exp1^(1/nu)*log(exp1))/nu^4 )
                   d2ldv2 <- -((1/nu^2) + (2*log(exp1))/nu^3 -d2exp2dv2)^2
                   d2ldv2 
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dd <- mu*(y^sigma)*log(y)
                   dexp1dm <- y^sigma
                   d2exp1dmdd <- (y^sigma)*log(y)
                   d2exp2dmdd <- (1/nu)*( ((1/nu)-1)*dexp1dd*dexp1dm*(exp1^((1/nu)-2)) + d2exp1dmdd*(exp1^((1/nu)-1)))
                   d2ldmdd <- -( (((1/nu)-1)/exp1^2)*(exp1*d2exp1dmdd - dexp1dm*dexp1dd) - d2exp2dmdd)^2
                   d2ldmdd
                 }, 
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dm <- y^sigma
                   d2exp2dmdv <- -dexp1dm*((1/nu^2)*(exp1^((1/nu)-1)) + (1/nu^3)*log(exp1)*(exp1^((1/nu)-1)))
                   d2ldmdv <- -( -( dexp1dm/(exp1*(nu^2)) + d2exp2dmdv ) )^2
                   d2ldmdv
                 }, 
                 
                 d2ldddv = function(y, mu, sigma) {
                   exp1 <- 1 + mu*(y^sigma)
                   dexp1dd <- mu*(y^sigma)*log(y)
                   d2exp2dddv <- -((1/nu^2)*dexp1dd*(exp1^((1/nu)-1)))*(1+ log(exp1)/nu)
                   d2ldddv <- -( -dexp1dd/(exp1*(nu^2)) -d2exp2dddv )^2
                   d2ldddv 
                 }, 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dGPW(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pGPW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname GPW
dGPW<-function(x,mu,sigma,nu,log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik <- log(mu) + log(sigma) - log(nu) + (sigma-1)*log(x) +
    ((1/nu)-1)*log(1+mu*(x^sigma)) + 
    (1-(1+mu*(x^sigma))^(1/nu))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname GPW
pGPW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- 1- exp(1-(1 + mu*(q^sigma))^(1/nu))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf 
}

#' @export
#' @rdname GPW
qGPW <- function(p,mu,sigma,nu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  term <- 1-log(1-p)
  q <- (((term^nu)-1) /mu)^(1/sigma)
  q
}

#' @export
#' @rdname GPW
rGPW <- function(n,mu,sigma,nu){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qGPW(p, mu,sigma,nu)
  r
}

#' @export
#' @rdname GPW
# Hazard function
hGPW<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dGPW(x,mu,sigma,nu, log = FALSE)/pGPW(q=x,mu,sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h
}
