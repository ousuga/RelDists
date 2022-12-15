#' The Beta Generalized Exponentiated family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Beta Generalized Exponentiated family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dBGE}
#' 
#' @details 
#' The Beta Generalized Exponentiated distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{\nu \tau}{B(\mu, \sigma)} \exp(-\nu x)(1- \exp(-\nu x))^{\tau \mu - 1} (1 - (1- \exp(-\nu x))^\tau)^{\sigma -1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a BGE distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BGE.R 
#'
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{barreto2010beta}{RelDists}
#' 
#'@importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BGE <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log"){
  mstats <- checklink("mu.link", "Beta Generalized Exponentiated", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Beta Generalized Exponentiated",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Beta Generalized Exponentiated", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Beta Generalized Exponentiated", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("BGE", "Beta Generalized Exponentiated"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                 nopar=4, 
                 type="Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link    = as.character(substitute(nu.link)), 
                 tau.link   = as.character(substitute(tau.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun    = vstats$linkfun,
                 tau.linkfun   = tstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv    = vstats$linkinv,
                 tau.linkinv   = tstats$linkinv, 
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr    = vstats$mu.eta, 
                 tau.dr   = tstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   dldm  <- tau*log(exp2) - digamma(mu) + digamma(mu+sigma)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   dldd  <- log(1 - exp2^tau) - digamma(sigma) + digamma(mu+sigma)
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau){
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldv  <- 1/nu - y + (tau*mu - 1)*y*exp1/exp2 - exp3*
                     tau*y*exp1*exp2^(tau-1)
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldt  <- 1/tau + mu*log(exp2) - exp3*log(exp2)*exp2^tau
                   dldt
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   dldm  <- tau*log(exp2) - digamma(mu) + digamma(mu+sigma)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   dldm  <- tau*log(exp2) - digamma(mu) + digamma(mu+sigma)
                   dldd  <- log(1 - exp2^tau) - digamma(sigma) + digamma(mu+sigma)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldm  <- tau*log(exp2) - digamma(mu) + digamma(mu+sigma)
                   dldv  <- 1/nu - y + (tau*mu - 1)*y*exp1/exp2 - exp3*
                     tau*y*exp1*exp2^(tau-1)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldm  <- tau*log(exp2) - digamma(mu) + digamma(mu+sigma)
                   dldt  <- 1/tau + mu*log(exp2) - exp3*log(exp2)*exp2^tau
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   dldd  <- log(1 - exp2^tau) - digamma(sigma) + digamma(mu+sigma)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldd  <- log(1 - exp2^tau) - digamma(sigma) + digamma(mu+sigma)
                   dldv  <- 1/nu - y + (tau*mu - 1)*y*exp1/exp2 - exp3*
                     tau*y*exp1*exp2^(tau-1)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldd  <- log(1 - exp2^tau) - digamma(sigma) + digamma(mu+sigma)
                   dldt  <- 1/tau + mu*log(exp2) - exp3*log(exp2)*exp2^tau
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldv  <- 1/nu - y + (tau*mu - 1)*y*exp1/exp2 - exp3*
                     tau*y*exp1*exp2^(tau-1)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldv  <- 1/nu - y + (tau*mu - 1)*y*exp1/exp2 - exp3*
                     tau*y*exp1*exp2^(tau-1)
                   dldt  <- 1/tau + mu*log(exp2) - exp3*log(exp2)*exp2^tau
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   exp1  <- exp(-nu*y)
                   exp2  <- 1 - exp1
                   exp3  <- (sigma - 1)/(1 - exp2^tau)
                   dldt  <- 1/tau + mu*log(exp2) - exp3*log(exp2)*exp2^tau
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dBGE(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres      = expression(rqres(pfun="pBGE", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial    = expression(nu    <- rep(1, length(y))),
                 tau.initial   = expression(tau   <- rep(1, length(y))), 
                 
                 mu.valid    = function(mu)    all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid    = function(nu)    all(nu > 0), 
                 tau.valid   = function(tau)   all(tau > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}




