#' The Four-Parameter Exponentiated Generalized Gamma family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Four-Parameter Exponentiated Generalized Gamma distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dFPEGG}
#' 
#' @details 
#' Four-Parameter Exponentiated Generalized Gamma distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = \frac {\nu \sigma}{\mu \tau(\tau)} (\frac {x}{\mu})^{\sigma \tau -1} exp\{(-{\frac {x}{\mu})^{\sigma}})\} 
#' \{ \gamma~_1~ (\tau, (\frac {x}{\mu})^{\sigma}) \}^{\nu -1} ,}
#' 
#' for x > 0. 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rFPEGG(n=200, mu=0.1, sigma=0.8, nu=10, tau=1.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family='FPEGG',
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
#' mu <- exp(-0.8 + -3 * x1)
#' sigma <- exp(0.77 - 2 * x2)
#' nu <- 10
#' tau <- 1.5
#' x <- rFPEGG(n=n, mu, sigma, nu, tau)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=FPEGG,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' exp(coef(mod, what="tau"))
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#'
#' \insertRef{cordeiro2011}{RelDists}
#' 
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
FPEGG <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log") {
  mstats <- checklink("mu.link", "Four-Parameter Exponentiated Generalized Gamma", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Four-Parameter Exponentiated Generalized Gamma",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Four-Parameter Exponentiated Generalized Gamma", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Four-Parameter Exponentiated Generalized Gamma", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("FPEGG", "Four-Parameter Exponentiated Generalized Gamma"), 
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
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "mu", delta = 1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / sigma + tau * log(t) - log(t) * t^sigma
                   B <- (nu - 1) * exp(-(t)^sigma) * ((t^sigma)^tau) * log(t)
                   C <-  zipfR::Igamma(tau, t^sigma)
                   dldd <- A + B / C
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / nu + log(zipfR::Igamma(tau, t^sigma) / base::gamma(tau))
                   dldv <- A 
                   dldv
                 },
                 
                 dldt = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "tau", delta = 1e-04)
                   dldt = as.vector(attr(nd, "gradient"))
                   dldt
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "mu", delta = 1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "mu", delta = 1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   t <- (y / mu)
                   A <- 1 / sigma + tau * log(t) - log(t) * t^sigma
                   B <- (nu - 1) * exp(-(t)^sigma) * ((t^sigma)^tau) * log(t)
                   C <-  zipfR::Igamma(tau, t^sigma)
                   dldd <- A + B / C
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "mu", delta = 1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   t <- (y / mu)
                   A <- 1 / nu + log(zipfR::Igamma(tau, t^sigma) / base::gamma(tau))
                   dldv <- A 
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "mu", delta = 1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "tau", delta = 1e-04)
                   dldt = as.vector(attr(nd, "gradient"))
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / sigma + tau * log(t) - log(t) * t^sigma
                   B <- (nu - 1) * exp(-(t)^sigma) * ((t^sigma)^tau) * log(t)
                   C <-  zipfR::Igamma(tau, t^sigma)
                   dldd <- A + B / C
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / sigma + tau * log(t) - log(t) * t^sigma
                   B <- (nu - 1) * exp(-(t)^sigma) * ((t^sigma)^tau) * log(t)
                   C <-  zipfR::Igamma(tau, t^sigma)
                   dldd <- A + B / C
                   D <- 1 / nu + log(zipfR::Igamma(tau, t^sigma) / base::gamma(tau))
                   dldv <- D
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / sigma + tau * log(t) - log(t) * t^sigma
                   B <- (nu - 1) * exp(-(t)^sigma) * ((t^sigma)^tau) * log(t)
                   C <-  zipfR::Igamma(tau, t^sigma)
                   dldd <- A + B / C
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "tau", delta = 1e-04)
                   dldt = as.vector(attr(nd, "gradient"))
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / nu + log(zipfR::Igamma(tau, t^sigma) / base::gamma(tau))
                   dldv <- A 
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   t <- (y / mu)
                   A <- 1 / nu + log(zipfR::Igamma(tau, t^sigma) / base::gamma(tau))
                   dldv <- A 
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "tau", delta = 1e-04)
                   dldt = as.vector(attr(nd, "gradient"))
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   nd = gamlss::numeric.deriv(dFPEGG(y, mu, sigma, nu, tau,
                                                      log = TRUE), "tau", delta = 1e-04)
                   dldt = as.vector(attr(nd, "gradient"))
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dFPEGG(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres = expression(rqres(pfun="pFPEGG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
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
