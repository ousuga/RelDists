#' The Flexible Weibull Extension distribution
#' 
#' @description 
#' The Flexible Weibull Extension distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @details 
#' The Flexible Weibull extension with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = (\mu + \sigma/x^2) exp(\mu x - \sigma/x) exp(-exp(\mu x-\sigma/x))}
#' 
#' for x>0.
#' 
#' @examples  
#' plot(1:5)
#' 
FWE <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link", "Flexible Weibull Extension", substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "Flexible Weibull Extension", substitute(sigma.link), c("log", "identity"))
  
  structure(list(family = c("FEW", "Flexible Weibull Extension"),
                 parameters = list(mu=TRUE, sigma=TRUE), 
                 nopar = 2, 
                 type = "Continuous",
                 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta,
                 
                 dldm = function(y, mu, sigma) (1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y)),
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dldm = 1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y)
                   d2ldm2 = -dldm * dldm
                 },
                 
                 dldd = function(y, mu, sigma) (1/(mu*y^2+sigma)-1/y+exp(mu*y-sigma/y)/y),
                 
                 d2ldd2 = function(y,mu,sigma) {
                   dldd = 1/(mu*y^2+sigma) - 1/y + exp(mu*y-sigma/y)/y
                   d2ldd2 = -dldd * dldd
                 },
                 
                 d2ldmdd = function(y,mu,sigma) {
                   dldm = 1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y)
                   dldd = 1/(mu*y^2+sigma) - 1/y + exp(mu*y-sigma/y)/y
                   d2ldmdd = -dldm * dldd
                   d2ldmdd
                 },
                 
                 G.dev.incr  = function(y, mu, sigma, ...) -2*dFWE(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pFWE", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(0.5, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  all(y > 0)
  ),
  class = c("gamlss.family","family"))
}
