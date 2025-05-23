#' The four parameter Exponentiated Generalized Gamma family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The four parameter Exponentiated Generalized Gamma distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dEGG}
#' 
#' @details 
#' Four parameter Exponentiated Generalized Gamma distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = \frac{\nu \sigma}{\mu \Gamma(\tau)} \left(\frac{x}{\mu}\right)^{\sigma \tau -1} \exp\left\{ - \left( \frac{x}{\mu} \right)^\sigma \right\} \left\{ \gamma_1\left( \tau, \left( \frac{x}{\mu} \right)^\sigma \right) \right\}^{\nu-1} ,}
#' 
#' for x > 0. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a EGG distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_EGG.R 
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Cordeiro, G. M., Ortega, E. M., & Silva, G. O. (2011). 
#' The exponentiated generalized gamma distribution with 
#' application to lifetime data. Journal of statistical 
#' computation and simulation, 81(7), 827-842.
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats pgamma
#' @importFrom VGAM pgamma.deriv
#' @export
EGG <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log") {
  mstats <- checklink("mu.link", "Exponentiated Generalized Gamma", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Generalized Gamma",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Exponentiated Generalized Gamma", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Exponentiated Generalized Gamma", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("EGG", "Exponentiated Generalized Gamma"), 
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
                   p1 <- -1/mu + (sigma*tau-1) * (-1/mu) 
                   p2 <- - sigma * (y/mu)^(sigma) * (-1/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * sigma * (y/mu)^(sigma) * (-1/mu)
                   dldm <- p1 + p2 + p3
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   p1 <- 1/sigma + tau * (log(y)-log(mu))
                   p2 <- -(y/mu)^sigma * log(y/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * (y/mu)^sigma * log(y/mu)
                   dldd <- p1 + p2 + p3
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau) {
                   aux <- (y/mu)^sigma
                   dldv <- 1/nu + log(pgamma(aux, tau))
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   p1 <- - digamma(tau) + sigma * (log(y)-log(mu))
                   aux <- (y/mu)^sigma
                   p2 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,3] / pgamma(aux, tau)
                   dldt <- p1 + p2
                   dldt
                 },
                 
                 # Segundas derivadas ---------------------------------
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   p1 <- -1/mu + (sigma*tau-1) * (-1/mu) 
                   p2 <- - sigma * (y/mu)^(sigma) * (-1/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * sigma * (y/mu)^(sigma) * (-1/mu)
                   dldm <- p1 + p2 + p3
                   
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   
                   p1 <- -1/mu + (sigma*tau-1) * (-1/mu) 
                   p2 <- - sigma * (y/mu)^(sigma) * (-1/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * sigma * (y/mu)^(sigma) * (-1/mu)
                   dldm <- p1 + p2 + p3
                   
                   p1 <- 1/sigma + tau * (log(y)-log(mu))
                   p2 <- -(y/mu)^sigma * log(y/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * (y/mu)^sigma * log(y/mu)
                   dldd <- p1 + p2 + p3
                   
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   
                   p1 <- -1/mu + (sigma*tau-1) * (-1/mu) 
                   p2 <- - sigma * (y/mu)^(sigma) * (-1/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * sigma * (y/mu)^(sigma) * (-1/mu)
                   dldm <- p1 + p2 + p3
                   
                   aux <- (y/mu)^sigma
                   dldv <- 1/nu + log(pgamma(aux, tau))
                   
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   
                   p1 <- -1/mu + (sigma*tau-1) * (-1/mu) 
                   p2 <- - sigma * (y/mu)^(sigma) * (-1/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * sigma * (y/mu)^(sigma) * (-1/mu)
                   dldm <- p1 + p2 + p3
                   
                   p1 <- - digamma(tau) + sigma * (log(y)-log(mu))
                   aux <- (y/mu)^sigma
                   p2 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,3] / pgamma(aux, tau)
                   dldt <- p1 + p2
                   
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   p1 <- 1/sigma + tau * (log(y)-log(mu))
                   p2 <- -(y/mu)^sigma * log(y/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * (y/mu)^sigma * log(y/mu)
                   dldd <- p1 + p2 + p3
                   
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   
                   p1 <- 1/sigma + tau * (log(y)-log(mu))
                   p2 <- -(y/mu)^sigma * log(y/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * (y/mu)^sigma * log(y/mu)
                   dldd <- p1 + p2 + p3
                   
                   aux <- (y/mu)^sigma
                   dldv <- 1/nu + log(pgamma(aux, tau))
                   
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   
                   p1 <- 1/sigma + tau * (log(y)-log(mu))
                   p2 <- -(y/mu)^sigma * log(y/mu)
                   aux <- (y/mu)^sigma
                   p3 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,1] / pgamma(aux, tau)
                   p3 <- p3 * (y/mu)^sigma * log(y/mu)
                   dldd <- p1 + p2 + p3
                   
                   p1 <- - digamma(tau) + sigma * (log(y)-log(mu))
                   aux <- (y/mu)^sigma
                   p2 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,3] / pgamma(aux, tau)
                   dldt <- p1 + p2
                   
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   aux <- (y/mu)^sigma
                   dldv <- 1/nu + log(pgamma(aux, tau)) 
                   
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   
                   aux <- (y/mu)^sigma
                   dldv <- 1/nu + log(pgamma(aux, tau))
                   
                   p1 <- - digamma(tau) + sigma * (log(y)-log(mu))
                   aux <- (y/mu)^sigma
                   p2 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,3] / pgamma(aux, tau)
                   dldt <- p1 + p2
                   
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   p1 <- - digamma(tau) + sigma * (log(y)-log(mu))
                   aux <- (y/mu)^sigma
                   p2 <- (nu-1) * VGAM::pgamma.deriv(aux, tau)[,3] / pgamma(aux, tau)
                   dldt <- p1 + p2
                   
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dEGG(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres = expression(rqres(pfun="pEGG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))),
                 sigma.initial = expression(sigma <- rep(1, length(y))),
                 nu.initial    = expression(nu    <- rep(1, length(y))),
                 tau.initial   = expression(tau   <- rep(1, length(y))),
                 
                 mu.valid    = function(mu)    all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid    = function(nu)    all(nu > 0),
                 tau.valid   = function(tau)   all(tau > 0),
                 y.valid     = function(y)     all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
