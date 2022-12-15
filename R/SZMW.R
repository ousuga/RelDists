#' The Sarhan and Zaindin's Modified Weibull family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Sarhan and Zaindin's Modified Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dSZMW}
#' 
#' @details 
#' The Sarhan and Zaindin's Modified Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=(\mu + \sigma \nu x^(\nu - 1)) \exp(- \mu x - \sigma x^\nu),}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a SZMW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_SZMW.R
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{sarhan2009modified}{RelDists}
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
SZMW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Sarhan and Zaindin's Modified Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Sarhan and Zaindin's Modified Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Sarhan and Zaindin's Modified Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("SZMW", "Sarhan and Zaindin's Modified Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
                       mu.link = as.character(substitute(mu.link)), 
                    sigma.link = as.character(substitute(sigma.link)), 
                       nu.link = as.character(substitute(nu.link)),
                 
                    mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                    nu.linkfun = vstats$linkfun,
                 
                    mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                    nu.linkinv = vstats$linkinv,
                 
                         mu.dr = mstats$mu.eta, 
                      sigma.dr = dstats$mu.eta, 
                         nu.dr = vstats$mu.eta,
                 
                 dldm = function(y, mu, sigma, nu) {
                   A    <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldm <- A - y
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   A    <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldd <- A * nu * (y^(nu - 1)) - y^nu
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   A    <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldv <- A*(sigma * y^(nu - 1))*(1+ nu*log(y)) - sigma*y^nu*log(y)
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A      <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldm   <- A - y
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd  = function(y, mu, sigma, nu) {
                   A       <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldm    <- A - y
                   dldd    <- A * nu * (y^(nu - 1)) - y^nu
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A       <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldm    <- A - y
                   dldv    <- A*(sigma * y^(nu - 1))*(1+ nu*log(y)) - sigma*y^nu*log(y)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   A      <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldd   <- A * nu * (y^(nu - 1)) - y^nu
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   A       <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldd    <- A * nu * (y^(nu - 1)) - y^nu
                   dldv    <- A*(sigma * y^(nu - 1))*(1+ nu*log(y)) - sigma*y^nu*log(y)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   A      <- 1/(mu + (sigma * nu * y^(nu - 1)))
                   dldv   <- A*(sigma * y^(nu - 1))*(1+ nu*log(y)) - sigma*y^nu*log(y)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                    G.dev.incr = function(y, mu, sigma, nu, ...) -2*dSZMW(y, mu, sigma, nu, log=TRUE), 
                         rqres = expression(rqres(pfun="pSZMW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                    nu.initial = expression(nu    <- rep(1, length(y))),
                 
                      mu.valid = function(mu)    all(mu > 0), 
                   sigma.valid = function(sigma) all(sigma > 0), 
                      nu.valid = function(nu)    all(nu > 0),
                 
                       y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}