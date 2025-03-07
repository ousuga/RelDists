#' The exponentiated XLindley family
#' 
#' @author Manuel Gutierrez Tangarife, \email{mgutierrezta@unal.edu.co}
#' 
#' @description 
#' The function \code{EXL()} defines The exponentiated XLindley, a two parameter 
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @references
#' Alomair, A. M., Ahmed, M., Tariq, S., Ahsan-ul-Haq, M., & Talib, J. (2024). An exponentiated XLindley distribution with properties, inference and applications. Heliyon, 10(3).
#' 
#' @seealso \link{EXL}.
#' 
#' @details 
#' The exponentiated XLindley with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#'  \deqn{f(x) = \frac{\sigma\mu^2(2+\mu + x)\exp(-\mu x)}{(1+\mu)^2}\left[1-
#'  \left(1+\frac{\mu x}{(1 + \mu)^2}\right) \exp(-\mu x)\right] ^ {\sigma-1} }
#' for \eqn{x \geq 0}, \eqn{\mu \geq 0} and \eqn{\sigma \geq 0}.
#' 
#' Note: In this implementation we changed the original parameters \eqn{\delta} for \eqn{\mu}
#' and \eqn{\alpha} for \eqn{\sigma} we did it to implement this distribution
#' within gamlss framework.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a EXL distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_EXL.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EXL <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "Exponentiated XLindley", substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "Exponentiated XLindley", substitute(sigma.link), c("log", "identity"))
  
  structure(list(family = c("EXL", "Exponentiated XLindley"),
                 parameters = list(mu=TRUE, sigma=TRUE), 
                 nopar = 2, 
                 type = "Continuous",
                 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta,
                 
                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dEXL(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00000001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                   
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dEXL(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00000001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   dd   <- gamlss::numeric.deriv(dEXL(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00000001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },
                 
                 d2ldd2 = function(y,mu,sigma) {
                   dd   <- gamlss::numeric.deriv(dEXL(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00000001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2 
                 },
                 
                 d2ldmdd = function(y,mu,sigma) {
                   dm   <- gamlss::numeric.deriv(dEXL(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00000001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dEXL(y, mu, sigma,log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00000001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },
                 
                 G.dev.incr = function(y, mu, sigma, ...) -2*dEXL(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pEXL", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 #mu.initial    = expression(    mu <- rep(initValuesEXL(y)[1], length(y)) ),     
                 #sigma.initial = expression( sigma <- rep(initValuesEXL(y)[2], length(y)) ), 
                 mu.initial    = expression(    mu <- rep(1, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu >= 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 y.valid = function(y)  all(y >= 0)
  ),
  class = c("gamlss.family","family"))
}
