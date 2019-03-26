#' The Modified Weibull family
#' 
#' @description
#' #' The Modified Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dMW}
#' 
#' @details 
#' The Modified Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \mu (\sigma + \nu x) x^(\sigma - 1) \exp(\nu x) \exp(-\mu x^(\sigma) \exp(\nu x)),}
#' 
#' for x > 0. 
#' 
#' @examples 
#'# Example 1
#'# Generating some random values with
#'# known mu, sigma and nu
#'y <- rMW(n=100, mu = 2, sigma = 1.5, nu = 0.2)
#'
#'# Fitting the model
#'require(gamlss)
#'
#'mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family= 'MW',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#'
#'# Extracting the fitted values for mu, sigma and nu
#'# using the inverse link function
#'exp(coef(mod, what='mu'))
#'exp(coef(mod, what='sigma'))
#'exp(coef(mod, what='nu'))
#'
#'# Example 2
#'# Generating random values under some model
#'n     <- 200
#'x1    <- rpois(n, lambda=2)
#'x2    <- runif(n)
#'mu    <- exp(3 -1 * x1)
#'sigma <- exp(2 - 2 * x2)
#'nu    <- 0.2
#'x     <- rMW(n=n, mu, sigma, nu)
#'
#'mod <- gamlss(x~x1, mu.fo=~x1, sigma.fo=~x2, nu.fo=~1, family=MW,
#'              control=gamlss.control(n.cyc=5000, trace=FALSE))
#'
#'coef(mod, what="mu")
#'coef(mod, what="sigma")
#'coef(mod, what='nu')
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MW <- function (mu.link="log", sigma.link="log", nu.link="log") 
{
  mstats <- checklink("mu.link", "Modified Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Modified Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Modified Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("MW", "Modified Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=4, 
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
                 
                 dldm = function(y, mu, sigma, nu) {
                   dldm  <- (1/mu)-y^sigma * exp(nu*y)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1  <- mu*(y^sigma)*exp(nu*y)
                   dldd  <- (1/(sigma+nu*y))+log(y)-log(y)*exp1
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1  <- mu*(y^sigma)*exp(nu*y)
                   dldv  <- y*((1/(sigma+nu*y))+1-exp1)
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   dldm   <- (1/mu)-y^sigma * exp(nu*y)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   dldm    <- (1/mu)-y^sigma * exp(nu*y)
                   exp1    <- mu*(y^sigma)*exp(nu*y)
                   dldd    <- (1/(sigma+nu*y))+log(y)-log(y)*exp1
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   dldm    <- (1/mu)-y^sigma * exp(nu*y)
                   exp1    <- mu*(y^sigma)*exp(nu*y)
                   dldv    <- y*((1/(sigma+nu*y))+1-exp1)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   exp1   <- mu*(y^sigma)*exp(nu*y)
                   dldd   <- (1/(sigma+nu*y))+log(y)-log(y)*exp1
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv =function(y, mu, sigma, nu) {
                   exp1    <- mu*(y^sigma)*exp(nu*y)
                   dldd    <- (1/(sigma+nu*y))+log(y)-log(y)*exp1
                   dldv    <- y*((1/(sigma+nu*y))+1-exp1)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1   <- mu*(y^sigma)*exp(nu*y)
                   dldv   <- y*((1/(sigma+nu*y))+1-exp1)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                    G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMW(y, mu, sigma, nu, log=TRUE), 
                         rqres = expression(rqres(pfun="pMW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                    nu.initial = expression(nu    <- rep(1, length(y))),
                 
                   mu.valid    = function(mu)    all(mu > 0), 
                   sigma.valid = function(sigma) all(sigma > 0), 
                   nu.valid    = function(nu)    all(nu > 0),
                 
                       y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
