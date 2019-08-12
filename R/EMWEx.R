#' The Exponentiated Modifien Weibull Extension family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Exponentiated Modifien Weibull Extension family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dEMWEx}
#' 
#' @details 
#' The Beta-Weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)= \nu \sigma \tau (\frac{x}{\mu})^{\sigma-1} \exp((\frac{x}{\mu})^\sigma +
#' \nu \mu (1- \exp((\frac{x}{\mu})^\sigma))) 
#' (1 - \exp (\nu\mu (1- \exp((\frac{x}{\mu})^\sigma))))^{\tau-1} ,}
#' 
#' for \eqn{x > 0}, \eqn{\nu> 0}, \eqn{\mu > 0}, \eqn{\sigma> 0} and \eqn{\tau > 0}. 
#' 
#' @examples
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rEMWEx(n=100, mu = 1, sigma =1.21, nu=1, tau=2)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family=EMWEx,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma, nu and tau
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' exp(coef(mod, what='nu'))
#' exp(coef(mod, what='tau'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 200
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(0.75 - x1)
#' sigma <- exp(0.5 - x2)
#' nu <- 1
#' tau <- 2
#' x <- rEMWEx(n=n, mu, sigma, nu, tau)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=EMWEx,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' exp(coef(mod, what="tau"))
#' 
#' @references
#'\insertRef{almalki2014modifications}{RelDists}
#'
#'\insertRef{sarhan2013exponentiated}{RelDists}
#' 
#'@importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EMWEx <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log"){
  mstats <- checklink("mu.link", "Beta-Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Beta-Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Beta-Weibull", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Beta-Weibull", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("EMWEx", "Beta-Weibull"), 
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
                   nm   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   nd   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau){
                   nv   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   nt   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=1e-04)
                   dldt <- as.vector(attr(nt, "gradient"))
                   dldt
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   nm   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   nm   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nd   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   nm   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nv   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   nm   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(nm, "gradient"))
                   nt   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=1e-04)
                   dldt <- as.vector(attr(nt, "gradient"))
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu, tau) {
                   nd   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   nd   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   nv   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   nd   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(nd, "gradient"))
                   nt   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=1e-04)
                   dldt <- as.vector(attr(nt, "gradient"))
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   nv   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   nv   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=1e-04)
                   dldv <- as.vector(attr(nv, "gradient"))
                   nt   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=1e-04)
                   dldt <- as.vector(attr(nt, "gradient"))
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   nt   <- gamlss::numeric.deriv(dEMWEx(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=1e-04)
                   dldt <- as.vector(attr(nt, "gradient"))
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dEMWEx(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres      = expression(rqres(pfun="pEMWEx", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial    = expression(nu    <- rep(1, length(y))),
                 tau.initial   = expression(tau   <- rep(1, length(y))), 
                 
                 mu.valid    = function(mu)    all(mu >= 0), 
                 sigma.valid = function(sigma) all(sigma >= 0), 
                 nu.valid    = function(nu)    all(nu >= 0), 
                 tau.valid   = function(tau)   all(tau >= 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}




