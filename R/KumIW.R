#' The Kumaraswamy Inverse Weibull family
#' 
#' @author Freddy Hernandez, \email{fhernanb@@unal.edu.co}
#' 
#' @description 
#' The Kumaraswamy Inverse Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dKumIW}
#' 
#' @details 
#' The Kumaraswamy Inverse Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \mu \sigma \nu x^{-\sigma - 1} \exp{- \mu x^{-\sigma}} (1 - \exp{- \mu x^{-\sigma}})^{\nu - 1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' The KumIW distribution with \eqn{\nu=1} corresponds with the IW distribution.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a KumIW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_KumIW.R
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{shahbaz2012kumaraswamy}{RelDists}
#'
#'@importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
KumIW <- function (mu.link="log", sigma.link="log", nu.link="log"){
  mstats <- checklink("mu.link", "Kumaraswamy Inverse-Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Kumaraswamy Inverse-Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Kumaraswamy Inverse-Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("KumIW", "Kumaraswamy Inverse-Weibull"), 
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
                 
                 # First derivates
                 dldm = function(y, mu, sigma, nu) {
                   dldm <- 1/mu - y^(-sigma) + (nu-1) * y^(-sigma) / (exp(mu*y^-sigma)-1)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   temp <- -(nu-1)*mu*log(y)*y^(-sigma) / (exp(mu*y^(-sigma))-1)
                   dldd <- 1/sigma - log(y) + mu*log(y)*y^(-sigma) + temp
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   dldv <- 1/nu + log(1-exp(-mu*y^(-sigma)))
                   dldv
                 },
                 
                 # Second derivates
                 d2ldm2 = function(y, mu, sigma, nu) {
                   dldm <- 1/mu - y^(-sigma) + (nu-1) * y^(-sigma) / (exp(mu*y^-sigma)-1)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   dldm <- 1/mu - y^(-sigma) + (nu-1) * y^(-sigma) / (exp(mu*y^-sigma)-1)
                   temp <- -(nu-1)*mu*log(y)*y^(-sigma) / (exp(mu*y^(-sigma))-1)
                   dldd <- 1/sigma - log(y) + mu*log(y)*y^(-sigma) + temp
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   dldm <- 1/mu - y^(-sigma) + (nu-1) * y^(-sigma) / (exp(mu*y^-sigma)-1)
                   dldv <- 1/nu + log(1-exp(-mu*y^(-sigma)))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   temp <- -(nu-1)*mu*log(y)*y^(-sigma) / (exp(mu*y^(-sigma))-1)
                   dldd <- 1/sigma - log(y) + mu*log(y)*y^(-sigma) + temp
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   temp <- -(nu-1)*mu*log(y)*y^(-sigma) / (exp(mu*y^(-sigma))-1)
                   dldd <- 1/sigma - log(y) + mu*log(y)*y^(-sigma) + temp
                   dldv <- 1/nu + log(1-exp(-mu*y^(-sigma)))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   dldv <- 1/nu + log(1-exp(-mu*y^(-sigma)))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dKumIW(y, mu, sigma, nu, log=TRUE), 
                 rqres      = expression(rqres(pfun="pKumIW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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




