#' The Marshall-Olkin Logistic-Exponential family
#' 
#' @description 
#' The Marshall-Olkin Logistic-Exponential family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dMOLE}
#' 
#' @details 
#' The Marshall-Olkin Logistic-Exponential distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#'
#' \eqn{f(x)= \frac{ \mu \sigma \nu \exp^{\nu x}[\exp^{\nu x}-1]^{-\mu-1}} {[1 + \sigma [\exp^{\nu x}-1]^{-\mu}]^2}}
#'
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 1}.
#'
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rMOLE(n=300, mu=0.6, sigma=5.5, nu=1)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='MOLE',
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
#' mu <- exp(-2.010 + 3 * x1)
#' sigma <- exp(2.704 -2 * x2)
#' nu <- 1
#' x <- rMOLE(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=MOLE,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' 
#' @references
#' \insertRef{mansoor2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MOLE <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Marshall-Olkin Logistic-Exponential", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall-Olkin Logistic-Exponential",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall-Olkin Logistic-Exponential", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("MOLE", "Marshall-Olkin Logistic-Exponential"), 
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
                   k <- exp(nu * y) - 1
                   B <- 1 / mu - log(k) + 2 * sigma * log(k) / (sigma + k^mu)    
                   dldm <- B
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   dldd <- 1 / sigma - 2 * k^(-mu) / (1 + sigma * k^(-mu))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   A <- 1 / nu + y + (-mu-1) * exp(nu * y) * y / k
                   dldv <- A + 2 * (sigma * k^(-mu-1) * mu *exp(nu * y) * y) / (1 + sigma * k^(-mu))
                   dldv
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   B <- 1 / mu - log(k) + 2 * sigma * log(k) / (sigma + k^mu)    
                   dldm <- B
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   B <- 1 / mu - log(k) + 2 * sigma * log(k) / (sigma + k^mu)    
                   dldm <- B
                   dldd <- 1 / sigma - 2 * k^(-mu) / (1 + sigma * k^(-mu))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   B <- 1 / mu - log(k) + 2 * sigma * log(k) / (sigma + k^mu)    
                   dldm <- B
                   A <- 1 / nu + y + (-mu-1) * exp(nu * y) * y / k
                   dldv <- A + 2 * (sigma * k^(-mu-1) * mu *exp(nu * y) * y) / (1 + sigma * k^(-mu))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   dldd <- 1 / sigma - 2 * k^(-mu) / (1 + sigma * k^(-mu))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   dldd <- 1 / sigma - 2 * k^(-mu) / (1 + sigma * k^(-mu))
                   A <- 1 / nu + y + (-mu-1) * exp(nu * y) * y / k
                   dldv <- A + 2 * (sigma * k^(-mu-1) * mu *exp(nu * y) * y) / (1 + sigma * k^(-mu))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   k <- exp(nu * y) - 1
                   A <- 1 / nu + y + (-mu-1) * exp(nu * y) * y / k
                   dldv <- A + 2 * (sigma * k^(-mu-1) * mu *exp(nu * y) * y) / (1 + sigma * k^(-mu))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMOLE(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pMOLE", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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
