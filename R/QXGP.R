#' The Quasi XGamma Poisson family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Quasi XGamma Poisson family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dQXGP}
#' 
#' @details 
#' The Quasi XGamma Poisson distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= K(\mu, \sigma, \nu)(\frac {\sigma^{2} x^{2}}{2} + \mu)
#'  exp(\frac{\nu exp(-\sigma x)(1 + \mu + \sigma x + \frac {\sigma^{2}x^{2}}{2})}{1+\mu} - \sigma x),}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 1}.
#' 
#' where
#' 
#' \eqn{K(\mu, \sigma, \nu) = \frac{\nu \sigma}{(exp(\nu)-1)(1+\mu)}} 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rQXGP(n=200, mu=4, sigma=2, nu=3)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='QXGP',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma and nu
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' exp(coef(mod, what='nu'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 2000
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(-2.19 + 3 * x1)
#' sigma <- exp(1 - 2 * x2)
#' nu <- 1
#' x <- rQXGP(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=QXGP,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' 
#' @references
#' \insertRef{subhradev2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
QXGP <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Quasi XGamma Poissonn", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Quasi XGamma Poisson",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Quasi XGamma Poisson", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("QXGP", "Quasi XGamma Poisson"), 
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
                   A    <- - 1 / (1 + mu) + 1 / ((1/2) * sigma * y^2 + mu)
                   B    <- nu * exp(-sigma * y) * (-sigma * y - (sigma^2 * y^2 / 2))
                   dldm <- A + B / (1 + mu)^2
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   C    <- 1 / sigma + sigma * y^2 / ((1/2) * sigma^2 * y^2 + mu)
                   D    <- (y + sigma * y^2) - y * (1 + mu + sigma * y + sigma^2 * y^2 / 2) 
                   dldd <- C + (D * nu * exp(-sigma * y) / (1 + mu)) - y
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   E    <- exp(-sigma * y) *(1 + mu + sigma * y + sigma^2 * y^2 / 2)
                   dldv <-  1 / nu - exp(nu) / (exp(nu) - 1) + E / (1 + mu)
                   dldv
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A    <- - 1 / (1 + mu) + 1 / ((1/2) * sigma * y^2 + mu)
                   B    <- nu * exp(-sigma * y) * (-sigma * y - (sigma^2 * y^2 / 2))
                   dldm <- A + B / (1 + mu)^2
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A    <- - 1 / (1 + mu) + 1 / ((1/2) * sigma * y^2 + mu)
                   B    <- nu * exp(-sigma * y) * (-sigma * y - (sigma^2 * y^2 / 2))
                   dldm <- A + B / (1 + mu)^2
                   C    <- 1 / sigma + sigma * y^2 / ((1/2) * sigma^2 * y^2 + mu)
                   D    <- (y + sigma * y^2) - y * (1 + mu + sigma * y + sigma^2 * y^2 / 2) 
                   dldd <- C + (D * nu * exp(-sigma * y) / (1 + mu)) - y
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A    <- - 1 / (1 + mu) + 1 / ((1/2) * sigma * y^2 + mu)
                   B    <- nu * exp(-sigma * y) * (-sigma * y - (sigma^2 * y^2 / 2))
                   dldm <- A + B / (1 + mu)^2
                   E    <- exp(-sigma * y) *(1 + mu + sigma * y + sigma^2 * y^2 / 2)
                   dldv <-  1 / nu - exp(nu) / (exp(nu) - 1) + E / (1 + mu)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   C    <- 1 / sigma + sigma * y^2 / ((1/2) * sigma^2 * y^2 + mu)
                   D    <- (y + sigma * y^2) - y * (1 + mu + sigma * y + sigma^2 * y^2 / 2) 
                   dldd <- C + (D * nu * exp(-sigma * y) / (1 + mu)) - y
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   C    <- 1 / sigma + sigma * y^2 / ((1/2) * sigma^2 * y^2 + mu)
                   D    <- (y + sigma * y^2) - y * (1 + mu + sigma * y + sigma^2 * y^2 / 2) 
                   dldd <- C + (D * nu * exp(-sigma * y) / (1 + mu)) - y
                   E    <- exp(-sigma * y) *(1 + mu + sigma * y + sigma^2 * y^2 / 2)
                   dldv <-  1 / nu - exp(nu) / (exp(nu) - 1) + E / (1 + mu)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   E    <- exp(-sigma * y) *(1 + mu + sigma * y + sigma^2 * y^2 / 2)
                   dldv <-  1 / nu - exp(nu) / (exp(nu) - 1) + E / (1 + mu)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dQXGP(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pQXGP", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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
