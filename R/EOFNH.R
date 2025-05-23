#' The Extended Odd Frechet-Nadarajah-Haghighi family
#' 
#' @author Helber Santiago Padilla, \email{hspadillar@unal.edu.co}
#' 
#' @description 
#' The Extended Odd Frechet-Nadarjad-Hanhighi family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @references
#' Nasiru, S. (2018). Extended Odd Fréchet‐G Family of Distributions 
#' Journal of Probability and Statistics, 2018(1), 2931326.
#' 
#' @seealso \link{dEOFNH}
#' 
#' @details 
#' The Extended Odd Frechet-Nadarajah-Haghighi distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \frac{\mu\sigma\nu\tau(1+\nu x)^{\sigma-1}e^{(1-(1+\nu x)^\sigma)}[1-(1-e^{(1-(1+\nu x)^\sigma)})^{\mu}]^{\tau-1}}{(1-e^{(1-(1+\nu x)^{\sigma})})^{\mu\tau+1}} e^{-[(1-e^{(1-(1+\nu x)^\sigma)})^{-\mu}-1]^{\tau}},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0}, \eqn{\nu > 0} and \eqn{\tau > 0}.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a EOFNH distribution in the \code{gamlss()} function.
#'  
#' @example examples/examples_EOFNH.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EOFNH <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log"){
  mstats <- checklink("mu.link", "Extended Odd Frechet-Nadarjad-Hanhighi", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Extended Odd Frechet-Nadarjad-Hanhighi",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Extended Odd Frechet-Nadarjad-Hanhighi", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Extended Odd Frechet-Nadarjad-Hanhighi", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("EOFNH", "Extended Odd Frechet-Nadarjad-Hanhighi"), 
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
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldm <- 1/mu - (tau-1)*((log(exp2)*exp2^mu)/exp3) - tau*log(exp2) 
                   #        + tau*log(exp2)*exp4^(tau-1)* exp2^(-mu)
                   nm   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldd <- 1/sigma + log(exp1) - exp1^sigma*log(exp1) 
                   #        - (tau-1)*((mu*exp(1-exp1^sigma)*log(exp1)*exp1^sigma*exp2^(mu-1))/exp3)
                   #        + tau*mu*exp(exp2)*log(exp1)*exp1^sigma*exp4^(tau-1)*exp2^(-mu-1) 
                   #        - ((sigma*y*exp(exp(1-exp1^sigma))*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nd   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau){
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldv <- 1/nu + (sigma-1)*(y/exp1) - y*sigma*exp1^(sigma-1) 
                   #        - (tau-1)*((y*sigma*mu*exp(exp2)*exp1^(sigma-1)*exp2^(mu-1))/exp3)
                   #        + sigma*mu*tau*y*exp(exp2)*exp4^(tau-1)*exp1^(sigma-1)*exp2^(-mu-1)
                   #        - ((sigma*y*exp(exp2)*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nv   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (1+nu*y)
                   exp2 <- 1- exp(1-exp1^sigma)
                   exp3 <- 1-exp2^mu
                   exp4 <- (exp2^(-mu)-1)
                   dldt <- 1/tau + log(exp3) - exp4^tau*log(exp4) - mu*log(exp2)
                   dldt
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldm <- 1/mu - (tau-1)*((log(exp2)*exp2^mu)/exp3) - tau*log(exp2) 
                   #        + tau*log(exp2)*exp4^(tau-1)* exp2^(-mu)
                   nm   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldm <- 1/mu - (tau-1)*((log(exp2)*exp2^mu)/exp3) - tau*log(exp2) 
                   #        + tau*log(exp2)*exp4^(tau-1)* exp2^(-mu)
                   #dldd <- 1/sigma + log(exp1) - exp1^sigma*log(exp1) 
                   #        - (tau-1)*((mu*exp(1-exp1^sigma)*log(exp1)*exp1^sigma*exp2^(mu-1))/exp3)
                   #        + tau*mu*exp(exp2)*log(exp1)*exp1^sigma*exp4^(tau-1)*exp2^(-mu-1) 
                   #        - ((sigma*y*exp(exp(1-exp1^sigma))*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nm   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nd   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldm <- 1/mu - (tau-1)*((log(exp2)*exp2^mu)/exp3) - tau*log(exp2) 
                   #        + tau*log(exp2)*exp4^(tau-1)* exp2^(-mu)
                   #dldv <- 1/nu + (sigma-1)*(y/exp1) - y*sigma*exp1^(sigma-1) 
                   #        - (tau-1)*((y*sigma*mu*exp(exp2)*exp1^(sigma-1)*exp2^(mu-1))/exp3)
                   #        + sigma*mu*tau*y*exp(exp2)*exp4^(tau-1)*exp1^(sigma-1)*exp2^(-mu-1)
                   #        - ((sigma*y*exp(exp2)*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nm   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nv   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (1+nu*y)
                   exp2 <- 1- exp(1-exp1^sigma)
                   exp3 <- 1-exp2^mu
                   exp4 <- (exp2^(-mu)-1)
                   #dldm <- 1/mu - (tau-1)*((log(exp2)*exp2^mu)/exp3) - tau*log(exp2) 
                   #        + tau*log(exp2)*exp4^(tau-1)* exp2^(-mu)
                   nm   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldt <- 1/tau + log(exp3) - exp4^tau*log(exp4) - mu*log(exp2)
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldd <- 1/sigma + log(exp1) - exp1^sigma*log(exp1) 
                   #        - (tau-1)*((mu*exp(1-exp1^sigma)*log(exp1)*exp1^sigma*exp2^(mu-1))/exp3)
                   #        + tau*mu*exp(exp2)*log(exp1)*exp1^sigma*exp4^(tau-1)*exp2^(-mu-1) 
                   #        - ((sigma*y*exp(exp(1-exp1^sigma))*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nd   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldd <- 1/sigma + log(exp1) - exp1^sigma*log(exp1) 
                   #        - (tau-1)*((mu*exp(1-exp1^sigma)*log(exp1)*exp1^sigma*exp2^(mu-1))/exp3)
                   #        + tau*mu*exp(exp2)*log(exp1)*exp1^sigma*exp4^(tau-1)*exp2^(-mu-1) 
                   #        - ((sigma*y*exp(exp(1-exp1^sigma))*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   #dldv <- 1/nu + (sigma-1)*(y/exp1) - y*sigma*exp1^(sigma-1) 
                   #        - (tau-1)*((y*sigma*mu*exp(exp2)*exp1^(sigma-1)*exp2^(mu-1))/exp3)
                   #        + sigma*mu*tau*y*exp(exp2)*exp4^(tau-1)*exp1^(sigma-1)*exp2^(-mu-1)
                   #        - ((sigma*y*exp(exp2)*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nd   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   nv   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (1+nu*y)
                   exp2 <- 1- exp(1-exp1^sigma)
                   exp3 <- 1-exp2^mu
                   exp4 <- (exp2^(-mu)-1)
                   ##dldd <- 1/sigma + log(exp1) - exp1^sigma*log(exp1) 
                   ##        - (tau-1)*((mu*exp(1-exp1^sigma)*log(exp1)*exp1^sigma*exp2^(mu-1))/exp3)
                   ##        + tau*mu*exp(exp2)*log(exp1)*exp1^sigma*exp4^(tau-1)*exp2^(-mu-1) 
                   ##        - ((sigma*y*exp(exp(1-exp1^sigma))*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nd   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   dldt <- 1/tau + log(exp3) - exp4^tau*log(exp4) - mu*log(exp2)
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (1+nu*y)
                   #exp2 <- 1- exp(1-exp1^sigma)
                   #exp3 <- 1-exp2^mu
                   #exp4 <- (exp2^(-mu)-1)
                   #dldv <- 1/nu + (sigma-1)*(y/exp1) - y*sigma*exp1^(sigma-1) 
                   #        - (tau-1)*((y*sigma*mu*exp(exp2)*exp1^(sigma-1)*exp2^(mu-1))/exp3)
                   #        + sigma*mu*tau*y*exp(exp2)*exp4^(tau-1)*exp1^(sigma-1)*exp2^(-mu-1)
                   #        - ((sigma*y*exp(exp2)*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nv   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (1+nu*y)
                   exp2 <- 1- exp(1-exp1^sigma)
                   exp3 <- 1-exp2^mu
                   exp4 <- (exp2^(-mu)-1)
                   #dldv <- 1/nu + (sigma-1)*(y/exp1) - y*sigma*exp1^(sigma-1) 
                   #        - (tau-1)*((y*sigma*mu*exp(exp2)*exp1^(sigma-1)*exp2^(mu-1))/exp3)
                   #        + sigma*mu*tau*y*exp(exp2)*exp4^(tau-1)*exp1^(sigma-1)*exp2^(-mu-1)
                   #        - ((sigma*y*exp(exp2)*exp1^(sigma-1)*(mu*tau+1))/exp2)
                   nv   <- gamlss::numeric.deriv(dEOFNH(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   dldt <- 1/tau + log(exp3) - exp4^tau*log(exp4) - mu*log(exp2)
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   exp1 <- (1+nu*y)
                   exp2 <- 1- exp(1-exp1^sigma)
                   exp3 <- 1-exp2^mu
                   exp4 <- (exp2^(-mu)-1)
                   dldt <- 1/tau + log(exp3) - exp4^tau*log(exp4) - mu*log(exp2)
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dEOFNH(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres      = expression(rqres(pfun="pEOFNH", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
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




