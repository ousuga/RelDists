#' The Marshall-Olkin Kappa family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Marshall-Olkin Kappa family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dMOK}
#' 
#' @details 
#' The Marshall-Olkin Kappa distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)=\frac{\tau\frac{\mu\nu}{\sigma}\left(\frac{x}{\sigma}\right)^{\nu-1} \left(\mu+\left(\frac{x}{\sigma}\right)^{\mu\nu}\right)^{-\frac{\mu+1}{\mu}}}{\left(\tau+(1-\tau)\left(\frac{\left(\frac{x}{\sigma}\right)^{\mu\nu}}{\mu+\left(\frac{x}{\sigma}\right)^{\mu\nu}}\right)^{\frac{1}{\mu}}\right)^2}}
#' 
#' for x > 0.
#' 
#' @example examples/examples_MOK.R
#'
#' @references
#'\insertRef{javed2018marshall}{RelDists}
#' 
#'@importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MOK <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log"){
  mstats <- checklink("mu.link", "Marshall-Olkin Kappa", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall-Olkin Kappa",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall-Olkin Kappa", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Marshall-Olkin Kappa", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("MOK", "Marshall-Olkin Kappa"), 
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
                   #exp1 <- (y/sigma)^(mu*nu)
                   #exp2 <- (exp1/(mu+exp1))^(1/mu)
                   #exp3 <- tau + (1-tau)*exp2
                   #dldm <- 1/mu + log(mu + exp1)/mu^2 - (nu*y^(mu*nu)*log(y/sigma)*(mu+1))/(mu*sigma^(mu*nu)*(mu+exp1)) - (2*(1-tau)*exp2*
                   #        ((((y/sigma)^(-mu*nu)*(exp1+mu)*((((nu*log(y/sigma)*exp1)/(mu+exp1))-((exp1*(nu*log(y/sigma)*exp1+ 1)))/(exp1+mu)^2)))/mu) -
                   #        (log(exp2))/mu^2))/exp3
                   nm   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (y/sigma)^(mu*nu)
                   #exp2 <- (exp1/(mu+exp1))^(1/mu)
                   #exp3 <- tau + (1-tau)*exp2
                   #dldd <- -nu/sigma + ((nu*y^mu*(mu+1))/(sigma^(mu*nu+1)*(mu+exp1))) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                   #         (((mu*nu*y*(y/sigma)^(2*mu*nu-1))/(sigma^2*(mu+exp1)^2))-((mu*nu*y*(y/sigma)^(mu*nu-1))/(sigma^2*(mu+exp1)))))/(mu*exp3))
                   nd   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau){
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   dldv <- 1/nu + log(y) - log(sigma) - ((exp1*log(y/sigma)*(mu+1))/(mu+exp1)) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                                                                                                     (((mu*log(y/sigma)*exp1)/(mu+exp1))-((mu*log(y/sigma)*(y/sigma)^(2*mu*nu))/(mu+exp1)^2)))/(mu*exp3))
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   dldt <- 1/tau - ((2*(1-exp2))/exp3)
                   dldt
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (y/sigma)^(mu*nu)
                   #exp2 <- (exp1/(mu+exp1))^(1/mu)
                   #exp3 <- tau + (1-tau)*exp2
                   #dldm <- 1/mu + log(mu + exp1)/mu^2 - (nu*y^(mu*nu)*log(y/sigma)*(mu+1))/(mu*sigma^(mu*nu)*(mu+exp1)) - (2*(1-tau)*exp2*
                   #       ((((y/sigma)^(-mu*nu)*(exp1+mu)*((((nu*log(y/sigma)*exp1)/(mu+exp1))-((exp1*(nu*log(y/sigma)*exp1+ 1)))/(exp1+mu)^2)))/mu) -
                   #       (log(exp2))/mu^2))/exp3 
                   nm   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (y/sigma)^(mu*nu)
                   #exp2 <- (exp1/(mu+exp1))^(1/mu)
                   #exp3 <- tau + (1-tau)*exp2
                   #dldm <- 1/mu + log(mu + exp1)/mu^2 - (nu*y^(mu*nu)*log(y/sigma)*(mu+1))/(mu*sigma^(mu*nu)*(mu+exp1)) - (2*(1-tau)*exp2*
                   #        ((((y/sigma)^(-mu*nu)*(exp1+mu)*((((nu*log(y/sigma)*exp1)/(mu+exp1))-((exp1*(nu*log(y/sigma)*exp1+ 1)))/(exp1+mu)^2)))/mu) -
                   #        (log(exp2))/mu^2))/exp3
                   #dldd <- -nu/sigma + ((nu*y^mu*(mu+1))/(sigma^(mu*nu+1)*(mu+exp1))) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                   #        (((mu*nu*y*(y/sigma)^(2*mu*nu-1))/(sigma^2*(mu+exp1)^2))-((mu*nu*y*(y/sigma)^(mu*nu-1))/(sigma^2*(mu+exp1)))))/(mu*exp3))
                   nm   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nd   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   #dldm <- 1/mu + log(mu + exp1)/mu^2 - (nu*y^(mu*nu)*log(y/sigma)*(mu+1))/(mu*sigma^(mu*nu)*(mu+exp1)) - (2*(1-tau)*exp2*
                   #        ((((y/sigma)^(-mu*nu)*(exp1+mu)*((((nu*log(y/sigma)*exp1)/(mu+exp1))-((exp1*(nu*log(y/sigma)*exp1+ 1)))/(exp1+mu)^2)))/mu) -
                   #        (log(exp2))/mu^2))/exp3 
                   nm   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldv <- 1/nu + log(y) - log(sigma) - ((exp1*log(y/sigma)*(mu+1))/(mu+exp1)) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                                                                                                     (((mu*log(y/sigma)*exp1)/(mu+exp1))-((mu*log(y/sigma)*(y/sigma)^(2*mu*nu))/(mu+exp1)^2)))/(mu*exp3))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   #dldm <- 1/mu + log(mu + exp1)/mu^2 - (nu*y^(mu*nu)*log(y/sigma)*(mu+1))/(mu*sigma^(mu*nu)*(mu+exp1)) - (2*(1-tau)*exp2*
                   #        ((((y/sigma)^(-mu*nu)*(exp1+mu)*((((nu*log(y/sigma)*exp1)/(mu+exp1))-((exp1*(nu*log(y/sigma)*exp1+ 1)))/(exp1+mu)^2)))/mu) -
                   #        (log(exp2))/mu^2))/exp3 
                   nm   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldt <- 1/tau - ((2*(1-exp2))/exp3)
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu, tau) {
                   #exp1 <- (y/sigma)^(mu*nu)
                   #exp2 <- (exp1/(mu+exp1))^(1/mu)
                   #exp3 <- tau + (1-tau)*exp2
                   #dldd <- -nu/sigma + ((nu*y^mu*(mu+1))/(sigma^(mu*nu+1)*(mu+exp1))) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                   #        (((mu*nu*y*(y/sigma)^(2*mu*nu-1))/(sigma^2*(mu+exp1)^2))-((mu*nu*y*(y/sigma)^(mu*nu-1))/(sigma^2*(mu+exp1)))))/(mu*exp3))
                   nd   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   #dldd <- -nu/sigma + ((nu*y^mu*(mu+1))/(sigma^(mu*nu+1)*(mu+exp1))) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                   #        (((mu*nu*y*(y/sigma)^(2*mu*nu-1))/(sigma^2*(mu+exp1)^2))-((mu*nu*y*(y/sigma)^(mu*nu-1))/(sigma^2*(mu+exp1)))))/(mu*exp3))
                   dldv <- 1/nu + log(y) - log(sigma) - ((exp1*log(y/sigma)*(mu+1))/(mu+exp1)) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                                                                                                     (((mu*log(y/sigma)*exp1)/(mu+exp1))-((mu*log(y/sigma)*(y/sigma)^(2*mu*nu))/(mu+exp1)^2)))/(mu*exp3))
                   nd   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   #dldd <- -nu/sigma + ((nu*y^mu*(mu+1))/(sigma^(mu*nu+1)*(mu+exp1))) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                   #        (((mu*nu*y*(y/sigma)^(2*mu*nu-1))/(sigma^2*(mu+exp1)^2))-((mu*nu*y*(y/sigma)^(mu*nu-1))/(sigma^2*(mu+exp1)))))/(mu*exp3))
                   nd   <- gamlss::numeric.deriv(dMOK(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   dldt <- 1/tau - ((2*(1-exp2))/exp3)
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   dldv <- 1/nu + log(y) - log(sigma) - ((exp1*log(y/sigma)*(mu+1))/(mu+exp1)) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                                                                                                     (((mu*log(y/sigma)*exp1)/(mu+exp1))-((mu*log(y/sigma)*(y/sigma)^(2*mu*nu))/(mu+exp1)^2)))/(mu*exp3))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   dldv <- 1/nu + log(y) - log(sigma) - ((exp1*log(y/sigma)*(mu+1))/(mu+exp1)) - ((2*(1-tau)*(exp1/(mu+exp1))^(1/mu-1)*
                                                                                                     (((mu*log(y/sigma)*exp1)/(mu+exp1))-((mu*log(y/sigma)*(y/sigma)^(2*mu*nu))/(mu+exp1)^2)))/(mu*exp3))
                   dldt <- 1/tau - ((2*(1-exp2))/exp3)
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   exp1 <- (y/sigma)^(mu*nu)
                   exp2 <- (exp1/(mu+exp1))^(1/mu)
                   exp3 <- tau + (1-tau)*exp2
                   dldt <- 1/tau - ((2*(1-exp2))/exp3)
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dMOK(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres      = expression(rqres(pfun="pMOK", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
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
