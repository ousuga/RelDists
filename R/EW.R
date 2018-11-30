#' The Exponentiated Weibull family
#' 
#' @description 
#' The Exponentiated Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dEW}
#' 
#' @details 
#' The Exponentiated Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}
#' 
#' for x > 0. 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rEW(n=100, mu=2, sigma=1.5, nu=0.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='EW',
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
#' n <- 20
#' x1 <- rpois(n, lambda=2)
#' x2 <- runif(n)
#' mu <- exp(2 + -3 * x1)
#' sigma <- exp(3 - 2 * x2)
#' nu <- 2
#' x <- rEW(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=EW, 
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' coef(mod, what="nu")
#' 
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EW <- function (mu.link="log", sigma.link="log", nu.link="log") 
{
  mstats <- checklink("mu.link",    "Exponentiated Weibull", 
                      substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Exponentiated Weibull", 
                      substitute(nu.link),    c("log", "own"))
  
  structure(list(family=c("EW", "Exponentiated Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
                 mu.link   =as.character(substitute(mu.link)), 
                 sigma.link=as.character(substitute(sigma.link)), 
                 nu.link   =as.character(substitute(nu.link)), 
                 
                 mu.linkfun   =mstats$linkfun, 
                 sigma.linkfun=dstats$linkfun, 
                 nu.linkfun   =vstats$linkfun, 
                 
                 mu.linkinv   =mstats$linkinv, 
                 sigma.linkinv=dstats$linkinv, 
                 nu.linkinv   =vstats$linkinv, 
                 
                 mu.dr   =mstats$mu.eta, 
                 sigma.dr=dstats$mu.eta, 
                 nu.dr   =vstats$mu.eta, 
                 
                 # mu
                 dldm = function(y, mu, sigma, nu) {
                   a <- mu * y^sigma
                   b <- 1 - exp(-a)
                   dldm <- 1/mu - y^sigma + (nu-1) * y^sigma * exp(-a) / b
                   dldm
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   nd = gamlss::numeric.deriv(dEW(y, mu, sigma, nu, log=TRUE), "mu", delta=1e-04)
                   dldm = as.vector(attr(nd, "gradient"))
                   d2ldm2 = -dldm * dldm
                   d2ldm2 
                 }, 
                 
                 # sigma
                 dldd = function(y, mu, sigma, nu) {
                   a <- mu * y^sigma
                   b <- 1 - exp(-a)
                   dldd <- 1/sigma + log(y) - a * log(y) + (nu-1) * a * log(y) * exp(-a) / b
                   dldd
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   nd = gamlss::numeric.deriv(dEW(y, mu, sigma, nu, log=TRUE), "sigma", delta=1e-04)
                   dldd = as.vector(attr(nd, "gradient"))
                   d2ldd2 = -dldd * dldd
                   d2ldd2 
                 }, 
                 
                 # nu
                 dldv = function(y, mu, sigma, nu) {
                   a <- mu * y^sigma
                   b <- 1 - exp(-a)
                   dldv <- 1/nu + log(b) 
                   dldv 
                 }, 
                 
                 d2ldv2 = function(nu) -1/nu^2, 
                 
                 # cross
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2exp1dmdd <- (y^sigma)*log(y)
                   d2exp2dmdd <- exp(-exp1)*(-dexp1dd*dexp1dm + d2exp1dmdd)
                   d2ldmdd <- -((y^sigma)*log(y) + 
                                  ((nu-1)/exp2^2)*(exp2*d2exp2dmdd - dexp2dm*dexp2dd))^2
                   d2ldmdd
                 },  
                 
                 d2ldmdv = function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   d2ldmdv <- -(dexp2dm/exp2)^2
                   d2ldmdv
                 },   
                 
                 d2ldddv = function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2ldddv <- -(dexp2dd/exp2)^2
                   d2ldddv 
                 },   
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dEW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pEW", type="Continuous",
                                          y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression( mu <-  rep(1, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 nu.initial = expression( nu <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >  0), 
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}