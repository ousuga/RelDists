#' The Weibull Poisson family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Weibull Poisson family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dWP}
#' 
#' @details 
#' The Weibull Poisson distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu e^{-\nu}} {1-e^{-\nu}} x^{\mu-1} \exp({-\sigma x^{\mu}+\nu \exp({-\sigma} x^{\mu}) }),}
#' 
#' for \eqn{x > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a WP distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_WP.R 
#' 
#' @references
#' Lu, Wanbo, and Daimin Shi. "A new compounding life distribution: 
#' the Weibull–Poisson distribution." Journal of applied 
#' statistics 39.1 (2012): 21-38.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
WP <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Weibull Poisson", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Weibull Poisson",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Weibull Poisson", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("WP", "Weibull Poisson"), 
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
                 
                 # Primeras derivadas ---------------------------------
                 dldm = function(y, mu, sigma, nu) {
                   A    <- nu * sigma * y^mu * exp(- sigma * y^mu) * log(y)
                   dldm <- 1 / mu + log(y) - A - sigma * y^mu * log(y)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   dldd <- 1 / sigma - nu * exp(- sigma * y^mu) * y^mu - y^mu
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   A    <- exp(-nu) / (1 - exp(-nu))
                   dldv <- 1 / nu - 1 - A + exp(- sigma * y^mu)
                   dldv
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A      <- nu * sigma * y^mu * exp(- sigma * y^mu) * log(y)
                   dldm   <- 1 / mu + log(y) - A - sigma * y^mu * log(y)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A       <- nu * sigma * y^mu * exp(- sigma * y^mu) * log(y)
                   dldm    <- 1 / mu + log(y) - A - sigma * y^mu * log(y)
                   dldd    <- 1 / sigma - nu * exp(- sigma * y^mu) * y^mu - y^mu
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A       <- nu * sigma * y^mu * exp(- sigma * y^mu) * log(y)
                   dldm    <- 1 / mu + log(y) - A - sigma * y^mu * log(y)
                   B       <- exp(-nu) / (1 - exp(-nu))
                   dldv    <- 1 / nu - 1 - B + exp(- sigma * y^mu)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   dldd   <- 1 / sigma - nu * exp(- sigma * y^mu) * y^mu - y^mu
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   dldd    <- 1 / sigma - nu * exp(- sigma * y^mu) * y^mu - y^mu
                   B       <- exp(-nu) / (1 - exp(-nu))
                   dldv    <- 1 / nu - 1 - B + exp(- sigma * y^mu)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   A      <- exp(-nu) / (1 - exp(-nu))
                   dldv   <- 1 / nu - 1 - A + exp(- sigma * y^mu)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dWP(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pWP", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression(mu       <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu       <- rep(1, length(y))),
                 
                 mu.valid = function(mu)       all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu)       all(nu > 0),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
