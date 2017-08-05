#' @name EW
#' 
#' @title 
#' The Exponentiated Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation  and hazard function for the exponentiated weibull  distribution with
#' parameters \code{mu}, \code{sigma} and \code{nu}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu scale parameter.
#' @param sigma,nu shape parameters.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The Exponentiated Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}
#' 
#' for x > 0. 
#' 
#' @return 
#' \code{dEW} gives the density, \code{pEW} gives the distribution 
#' function, \code{qEW} gives the quantile function, \code{rEW}
#' generates random deviates and \code{hEW} gives the hazard function.
#'
#' @export
#' @examples  
#' ## The probability density function
#' curve(dEW(x, mu = 2, sigma = 1.5, nu = 0.5), from = 0, to = 2, ylim = c(0, 2.5), col = "red", las = 1, ylab = "The probability density function") 
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pEW(x, mu = 2, sigma = 1.5, nu = 0.5), from = 0, to = 2,  col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pEW(x, mu = 2, sigma = 1.5, nu = 0.5, lower.tail = FALSE), from = 0, to = 2,  col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qEW(p, mu = 2, sigma = 1.5, nu = 0.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pEW(x, mu = 2, sigma = 1.5, nu = 0.5),  from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rEW(n = 10000, mu = 2, sigma = 1.5, nu = 0.5), freq = FALSE, xlab = "x", las = 1, main = "")
#' curve(dEW(x, mu = 2, sigma = 1.5, nu = 0.5),  from = 0, add = TRUE, col = "red") 
#' 
#' ## The Hazard function
#' curve(hEW(x, mu = 2,sigma = 1.5, nu = 0.5), from = 0, to = 2, ylim = c(0, 7), col = "red", ylab = "The Hazard function")
#' 
EW <- function (mu.link = "log", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link",    "Exponentiated Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Exponentiated Weibull", substitute(nu.link),    c("log", "own"))
  
  structure(list(family = c("EW", "Exponentiated Weibull"), 
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
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- (y^sigma)
                   dexp2dm <- exp(-exp1)*dexp1dm
                   dldm <- 1/mu - y^sigma + ((nu-1)*dexp2dm)/exp2
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- (y^sigma)
                   dexp2dm <-  exp(-exp1)*dexp1dm
                   d2exp2dm2 <- -exp(-exp1)*(dexp1dm)^2
                   d2ldm2 <- -(-(1/mu^2) + ((nu-1)/exp2^2)*(exp2*d2exp2dm2-dexp2dm^2))^2
                   d2ldm2 
                 }, 
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   dldd <- 1/sigma + log(y) -exp1*log(y) + ((nu-1)*dexp2dd)/exp2
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1) 
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2exp1dd2 <- log(y)*dexp1dd   
                   d2exp2dd2 <- exp(-exp1)*(d2exp1dd2-dexp1dd^2)
                   d2ldd2 <- -(-(1/sigma^2) -log(y)*dexp1dd - ((nu-1)/exp2^2)*(exp2*d2exp2dd2-dexp2dd^2))^2
                   d2ldd2
                 }, 
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dldv <- 1/nu + log(exp2) 
                   dldv 
                 }, 
                 
                 d2ldv2 = function(nu) -(-(1/nu^2))^2, 
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2exp1dmdd <- (y^sigma)*log(y)
                   d2exp2dmdd <- exp(-exp1)*(-dexp1dd*dexp1dm + d2exp1dmdd)
                   d2ldmdd <- -((y^sigma)*log(y) + ((nu-1)/exp2^2)*(exp2*d2exp2dmdd - dexp2dm*dexp2dd))^2
                   d2ldmdd
                 },  
                 
                 d2ldmdv = function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   d2ldmdv <- -(dexp2dm/exp2)^2
                   d2ldmdv
                 },   
                 
                 d2ldddv = function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2ldddv <- -(dexp2dd/exp2)^2
                   d2ldddv 
                 },   
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dEW(y, mu, sigma, nu, log = TRUE), 
                 rqres = expression(rqres(pfun = "pEW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                 mu.initial = expression( mu <-  rep(1, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 nu.initial = expression( nu <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}
#' @export
#' @rdname EW
dEW<-function(x,mu,sigma,nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  loglik<- log(nu) +log(mu) + log(sigma) + (sigma-1)*log(x)-
    mu*(x^sigma) + (nu-1)*log(1-exp(-mu*(x^sigma)))
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname EW
pEW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  cdf <- (1-exp(-mu*(q^sigma)))^nu
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname EW
qEW <- function(p,mu,sigma,nu, lower.tail = TRUE, log.p = FALSE){
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
  
  q <- ((-1/mu)*log(1-p^(1/nu)))^(1/sigma)
  q
}
#' @export
#' @rdname EW
rEW <- function(n,mu,sigma, nu){
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
  r <- qEW(p, mu,sigma,nu)
  r
}
#' @export
#' @rdname EW
hEW<-function(x,mu,sigma, nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu<=0)) 
    stop(paste("nu must be positive", "\n", ""))
  
  h <- dEW(x,mu,sigma, nu, log = FALSE)/pEW(q=x,mu,sigma, nu, lower.tail=FALSE, log.p = FALSE)
  h
}

