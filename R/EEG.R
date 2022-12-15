#' The Extended Exponential Geometric family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Extended Exponential Geometric family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma. 
#' 
#' @seealso \link{dEEG}
#' 
#' @details 
#' The Extended Exponential Geometric distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x)= \mu \sigma \exp(-\mu x)(1 - (1 - \sigma)\exp(-\mu x))^{-2},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0} and \eqn{\sigma > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a EEG distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_EEG.R 
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#'
#' \insertRef{adamidis2005extension}{RelDists}
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EEG <- function (mu.link="log", sigma.link="log"){
  mstats <- checklink("mu.link", "Extended Exponential Geometric", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Extended Exponential Geometric",
                      substitute(sigma.link), c("log", "own"))
  
  structure(list(family=c("EEG", "Extended Exponential Geometric"), 
                 parameters=list(mu=TRUE, sigma=TRUE), 
                 nopar=2, 
                 type="Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta,
                 
                 dldm = function(y, mu, sigma) {
                   exp1 <- (1-sigma)*exp(-mu*y)
                   exp2 <- 1/1-exp1
                   dldm <- 1/mu - y - 2*y*exp1*exp2
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   exp1 <- (1-sigma)*exp(-mu*y)
                   exp2 <- 1/1-exp1
                   dldd <- 1/sigma - 2*exp(-mu*y)*exp2
                   dldd
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   exp1   <- (1-sigma)*exp(-mu*y)
                   exp2   <- 1/1-exp1
                   dldm   <- 1/mu - y - 2*y*exp1*exp2
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   exp1    <- (1-sigma)*exp(-mu*y)
                   exp2    <- 1/1-exp1
                   dldm    <- 1/mu - y - 2*y*exp1*exp2
                   dldd    <- 1/sigma - 2*exp(-mu*y)*exp2
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 
                 d2ldd2  = function(y, mu, sigma) {
                   exp1   <- (1-sigma)*exp(-mu*y)
                   exp2   <- 1/1-exp1
                   dldd   <- 1/sigma - 2*exp(-mu*y)*exp2
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, ...) -2*dEEG(y, mu, sigma, log=TRUE), 
                 rqres      = expression(rqres(pfun="pEEG", type="Continuous", y=y, mu=mu, sigma=sigma)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))),
                 
                 mu.valid    = function(mu)    all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}




