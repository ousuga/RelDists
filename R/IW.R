#' The Inverse Weibull family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Inverse Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @seealso \link{dIW}
#' 
#' @details 
#' The Inverse Weibull distribution with parameters \code{mu}, 
#' \code{sigma} has density given by
#' 
#' \eqn{f(x) = \mu \sigma x^{-\sigma-1} \exp(\mu x^{-\sigma})}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0} and \eqn{\sigma > 0} 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a IW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_IW.R 
#'
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Drapella, A. (1993). The complementary Weibull distribution: 
#' unknown or just forgotten?. Quality and reliability engineering 
#' international, 9(4), 383-385.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
IW <- function (mu.link="log", sigma.link="log"){
   mstats <- checklink("mu.link", "Inverse Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Inverse Weibull",
                      substitute(sigma.link), c("log", "own"))
  
  structure(list(family=c("IW", "Inverse Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE), 
                 nopar=2, 
                 type="Continuous", 
                 
                 mu.link       = as.character(substitute(mu.link)), 
                 sigma.link    = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 
                 mu.dr         = mstats$mu.eta, 
                 sigma.dr      = dstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma) {
                   dldm <- 1/mu - y^(-sigma)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dldd <- 1/sigma - log(y)+ mu*y^(-sigma)*log(y)
                   dldd
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dldm   <- 1/mu - y^(-sigma)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 
                 },
                 
                 d2ldd2 = function(y, mu, sigma) {
                   dldd   <- 1/sigma - log(y)+ mu*y^(-sigma)*log(y)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   dldm    <- 1/mu - y^(-sigma)
                   dldd    <- 1/sigma - log(y) +  mu*y^(-sigma)*log(y)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 G.dev.incr    = function(y, mu, sigma,...) -2*dIW(y, mu, sigma, log=TRUE), 
                 rqres         = expression(rqres(pfun="pIW", type="Continuous",
                                                  y=y, mu=mu, sigma=sigma)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 
                 mu.valid      = function(mu)    all(mu > 0), 
                 sigma.valid   = function(sigma) all(sigma > 0), 
                 
                 y.valid       = function(y)     all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}