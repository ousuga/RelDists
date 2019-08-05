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
#' @examples 
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rEEG(n=100, mu = 1, sigma =1.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, family=EEG,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma, nu and tau
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 200
#' x1 <- runif(n, min=0.1, max=0.2)
#' x2 <- runif(n, min=0.1, max=0.15)
#' mu <- exp(0.75 - x1)
#' sigma <- exp(0.5 - x2)
#' x <- rEEG(n=n, mu, sigma)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, family=EEG,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
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




