#' The Additive Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Additive Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dAddW}
#' 
#' @details 
#' Additive Weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = (\mu\nu x^{\nu - 1} + \sigma\tau x^{\tau - 1}) \exp({-\mu x^{\nu} - \sigma x^{\tau} }),}
#' 
#' for \eqn{x > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a 
#' AddW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_AddW.R 
#' 
#' @references
#' Almalki, S. J. (2018). A reduced new modified Weibull distribution. 
#' Communications in Statistics-Theory and Methods, 47(10), 2297-2313.
#'
#' Xie, M., & Lai, C. D. (1996). Reliability analysis using an additive 
#' Weibull model with bathtub-shaped failure rate function. 
#' Reliability Engineering & System Safety, 52(1), 87-93.
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
AddW <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log") {
  mstats <- checklink("mu.link", "Additive Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Additive Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Additive Weibull", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Additive Weibull", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("AddW", "Additive Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                 nopar=4, 
                 type="Continuous", 
                 
                 mu.link = as.character(substitute(mu.link)), 
              sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)),
                tau.link = as.character(substitute(tau.link)),
                 
              mu.linkfun = mstats$linkfun, 
           sigma.linkfun = dstats$linkfun, 
              nu.linkfun = vstats$linkfun,
             tau.linkfun = tstats$linkfun,
                 
              mu.linkinv = mstats$linkinv, 
           sigma.linkinv = dstats$linkinv, 
              nu.linkinv = vstats$linkinv,
             tau.linkinv = tstats$linkinv,
                 
                   mu.dr = mstats$mu.eta, 
                sigma.dr = dstats$mu.eta, 
                   nu.dr = vstats$mu.eta,
                  tau.dr = tstats$mu.eta,
                 
           # Primeras derivadas ---------------------------------
           dldm = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm <- (nu * y^(nu - 1))/A - y^nu
             dldm
           },
           
           dldd = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd <- (tau * y^(tau - 1))/A - y^tau
             dldd
           },
           
           dldv = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B    <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv <- (mu * B)/A - mu * log(y) * y^nu 
             dldv
           },
           
           dldt = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B    <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt <- (sigma * B)/A - sigma * log(y) * y^tau 
             dldt
           },
           
           # Segundas derivadas ---------------------------------
           d2ldm2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm   <- (nu * y^(nu - 1))/A - y^nu
             d2ldm2 <- -dldm * dldm
             d2ldm2
           },
           
           d2ldmdd = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             d2ldmdd <- -dldm * dldd
             d2ldmdd
           },
           
           d2ldmdv = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldmdv <- -dldm * dldv
             d2ldmdv
           },
           
           d2ldmdt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             B       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldmdt <- -dldm * dldt
             d2ldmdt
           },
           
           d2ldd2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd   <- (tau * y^(tau - 1))/A - y^tau
             d2ldd2 <- -dldd * dldd
             d2ldd2
           },
           
           d2ldddv = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldddv <- -dldd * dldv
             d2ldddv
           },
           
           d2ldddt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             B       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldddt <- -dldd * dldt
             d2ldddt
           },
           
           d2ldv2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B      <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv   <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldv2 <- -dldv * dldv
             d2ldv2
           },
           
           d2ldvdt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             C       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * C)/A - sigma * log(y) * y^tau 
             d2ldvdt <- -dldv * dldt
             d2ldvdt
           },
           
           d2ldt2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B      <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt   <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldt2 <- -dldt * dldt
             d2ldt2
           },
           
            G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dAddW(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres = expression(rqres(pfun="pAddW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
            mu.initial = expression(mu    <- rep(1, length(y))), 
         sigma.initial = expression(sigma <- rep(1, length(y))), 
            nu.initial = expression(nu    <- rep(1, length(y))),
           tau.initial = expression(tau    <- rep(1, length(y))),
                 
              mu.valid = function(mu)    all(mu > 0), 
           sigma.valid = function(sigma) all(sigma > 0), 
              nu.valid = function(nu)    all(nu > 0),
             tau.valid = function(tau)    all(tau > 0),
                 
               y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}