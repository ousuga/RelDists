#' The Gamma Weibull family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Gamma Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dGammaW}
#' 
#' @details 
#' The Gamma Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \frac{\sigma \mu^{\nu}}{\Gamma (\nu)} x^{\nu \sigma - 1} \exp(-\mu x^\sigma),}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a GammaW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_GammaW.R
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Stacy, E. W. (1962). A generalization of the gamma 
#' distribution. The Annals of mathematical statistics, 1187-1192.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GammaW <- function (mu.link="log", sigma.link="log", nu.link="log") 
{
  mstats <- checklink("mu.link", "Gamma Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Gamma Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Gamma Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("GammaW", "Gamma Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
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
                   dldm <- (nu/mu) - y^sigma
                   dldm
                 },
                 
                 dldd    =  function(y, mu, sigma, nu) {
                   dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                   dldd
                 },
                 
                 dldv    =  function(y, mu, sigma, nu){
                   dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                 dldv
                 },
                 
                d2ldm2   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  d2ldm2 <- -dldm * dldm
                  d2ldm2
                  },
                
                d2ldmdd   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  d2ldmdd <- -dldm * dldd
                  d2ldmdd
                  },
                
                d2ldmdv   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldmdv <- -dldm * dldv
                  d2ldmdv
                  },
                
                d2ldd2   =  function(y, mu, sigma, nu, tau) {
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  d2ldd2 <- -dldd * dldd
                  d2ldd2
                  },
                
                d2ldddv   =  function(y, mu, sigma, nu, tau) {
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldddv <- -dldd * dldv
                  d2ldddv
                  },
                
                d2ldv2   =  function(y, mu, sigma, nu, tau) {
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldv2 <- -dldv * dldv
                  d2ldv2
                  },
                
                
                G.dev.incr = function(y, mu, sigma, nu, ...) -2*dGammaW(y, mu, sigma, nu, log=TRUE), 
                rqres      = expression(rqres(pfun="pGammaW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                
                mu.initial    = expression(mu    <- rep(1, length(y))), 
                sigma.initial = expression(sigma <- rep(1, length(y))), 
                nu.initial    = expression(nu    <- rep(1, length(y))),
                
                mu.valid    = function(mu)    all(mu > 0), 
                sigma.valid = function(sigma) all(sigma > 0), 
                nu.valid    = function(nu)    all(nu > 0),
                
                y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}

