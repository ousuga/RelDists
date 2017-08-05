#' @name WG
#' 
#' 
#' @title 
#' The Weibull Geometric Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the weibull geometric distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu shape parameter.
#' @param sigma scale parameter.
#' @param nu parameter of geometric random variable.             
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'  
#' @details 
#' The weibull geometric distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#' 
#' f(x) = (mu*(sigma)^mu*(1-nu)*x^(mu-1)*exp(-(sigma*x)^mu))/(1- nu*exp(-(sigma*x)^mu))^2
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dWG} gives the density, \code{pWG} gives the distribution 
#' function, \code{qWG} gives the quantile function, \code{rWG}
#' generates random deviates and \code{hWG} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dWG(x, mu = 0.5, sigma = 0.2, nu = 0.95), from = 0, to = 5, ylim = c(0, 1.5), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pWG(x, mu = 0.5, sigma = 0.2, nu = 0.95), from = 0, to = 5, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pWG(x, mu = 0.5, sigma = 0.2, nu = 0.95, lower.tail = FALSE), from = 0, to = 5, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qWG(p = p, mu = 0.5, sigma = 0.2, nu = 0.95), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pWG(x, mu = 0.5, sigma = 0.2, nu = 0.95), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rWG(1000, mu = 0.5, sigma = 0.2, nu = 0.95), freq = FALSE, xlab = "x", ylim = c(0, 0.2), las = 1, main = "")
#' curve(dWG(x, mu = 0.5, sigma = 0.2, nu = 0.95),  from = 0, add = TRUE, col = "red", ylim = c(0, 0.2))
#' 
#' ## The Hazard function
#' curve(hWG(x, mu = 0.5, sigma = 0.2, nu = 0.95), from = 0, to = 15, ylim = c(0, 2.2), col = "red", ylab = "The hazard function", las = 1)
WG <- function (mu.link = "log", sigma.link = "log", nu.link = "logit") 
{
  mstats <- checklink("mu.link",    "Weibull Geometric", substitute(mu.link),    c("identity", "own"))
  dstats <- checklink("sigma.link", "Weibull Geometric", substitute(sigma.link), c("identity", "own"))
  vstats <- checklink("nu.link",    "Weibull Geometric", substitute(nu.link),    c("logit","probit", "own"))
  
  structure(list(family = c("WG", "Weibull Geometric"), 
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
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dm <- exp1*log(sigma*y)
                   dexp2dm <- -exp2*dexp1dm
                   dldm <- 1/mu + log(sigma) + log(y) -dexp1dm + (2*dexp2dm)/(1-exp2)
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dm <- exp1*log(sigma*y)
                   dexp2dm <- -exp2*dexp1dm
                   d2exp1dm2 <- exp1*(log(sigma*y))^2
                   d2exp2dm2 <- (exp1^2)*exp2*(log(sigma*y))^3
                   d2ldm2 <- -(-(1/mu^2) - d2exp1dm2 + (2/(1-exp2)^2)*((1-exp2)*d2exp2dm2 + dexp2dm^2))^2
                   d2ldm2
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dd <- mu*y*(sigma*y)^(mu-1)
                   dexp2dd <- -exp2*dexp1dd
                   dldd <- mu/sigma -dexp1dd  + (2*dexp2dd)/(1-exp2) 
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dd <- mu*y*(sigma*y)^(mu-1)
                   dexp2dd <- -exp2*dexp1dd
                   d2exp1dd2 <- mu*(mu-1)*(y^2)*(sigma*y)^(mu-2)
                   d2exp2dd2 <- exp2*dexp1dd*d2exp1dd2
                   d2ldd2 <- -(-(mu/sigma^2) -d2exp1dd2 +  (2/(1-exp2)^2)*((1-exp2)*d2exp2dd2+dexp2dd^2))^2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dldv <- -(1/(1-nu)) + 2/(1-exp2)
                   dldv 
                 }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp2dv <- exp(-exp1)
                   d2ldv2 <- -(-(1/(1-nu)^2) + (2*dexp2dv)/(1-exp2)^2 )^2
                   d2ldv2 
                 }, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dm <- exp1*log(sigma*y)
                   dexp2dm <- -exp2*dexp1dm
                   dexp1dd <- mu*y*(sigma*y)^(mu-1)
                   dexp2dd <- -exp2*dexp1dd
                   d2exp1dmdd <- log(sigma*y)*dexp1dd + exp1/sigma
                   d2exp2dmdd <- -dexp2dd*dexp1dm - exp2*d2exp1dmdd
                   d2ldmdd <- -(1/sigma -d2exp1dmdd + (2/(1-exp2)^2)*((1-exp2)*d2exp2dmdd + dexp2dm*dexp2dd))^2
                   d2ldmdd
                 },  
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dm <- exp1*log(sigma*y)
                   dexp2dm <- -exp2*dexp1dm
                   dexp2dv <- exp(-exp1)
                   d2ldmdv <- -(2*dexp2dm*dexp2dv/(1-exp2)^2)^2
                   d2ldmdv
                 },   
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- (sigma*y)^mu
                   exp2 <- nu*exp(-exp1)
                   dexp1dd <- mu*y*(sigma*y)^(mu-1)
                   dexp2dd <- -exp2*dexp1dd
                   dexp2dv <- exp(-exp1)
                   d2exp2dddv <- -dexp2dv*dexp1dd
                   d2ldddv <- -((2/(1-exp2)^2)*((1-exp2)*d2exp2dddv + dexp2dd*dexp2dv))^2
                   d2ldddv 
                 },   
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dWG(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pWG", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0 & nu < 1), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname WG

dWG<-function(x,mu,sigma,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <=0  | nu>=1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  loglik<- log(mu) + mu*log(sigma) + log(1-nu) + (mu-1)*log(x) - 
    (sigma*x)^mu - 2*log(1-nu*exp(-(sigma*x)^mu))
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname WG
pWG <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <=0  | nu>=1  )) 
    stop(paste("p must be between zero and one", "\n", ""))
  
  cdf <- (1-exp(-(sigma*q)^mu))/(1-nu*exp(-(sigma*q)^mu))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname WG
qWG <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <=0  | nu>=1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  
  fda <- function(x,mu, sigma,nu){
    (1- exp(-(sigma*x)^mu))/(1-(nu*exp(-(sigma*x)^mu)))
  }
  
  fda1 <- function(x, mu, sigma, nu, p) {fda(x, mu, sigma,nu) - p}
  
  r_de_la_funcion <- function(mu, sigma, nu,p) {
    uniroot(fda1, interval=c(0,1e+06), mu, sigma, nu,p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu,p)
  q
  
}

#' @export
#' @rdname WG
rWG <- function(n,mu,sigma,nu){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <=0  | nu>=1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qWG(p, mu,sigma,nu)
  r
}
#' @export
#' @rdname WG
hWG<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <=0  | nu>=1  )) 
    stop(paste("nu must be between zero and one", "\n", ""))
  
  h <- dWG(x,mu,sigma,nu, log = FALSE)/pWG(q=x,mu,sigma,nu, lower.tail=FALSE, log.p = FALSE)
  h  
}
