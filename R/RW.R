#' The Reflected Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' Reflected Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @seealso \link{dRW}
#' 
#' @details 
#' The Reflected Weibull Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(y) = \mu\sigma (-y) ^{\sigma - 1} e ^ {-\mu(-y)^\sigma},}
#' 
#' for y < 0
#' 
#' @returns Returns a gamlss.family object which can be used to fit a RW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_RW.R
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#' 
#' Cohen, A. C. (1973). The reflected Weibull distribution. 
#' Technometrics, 15(4), 867-873.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
RW <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "Reflected Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Reflected Weibull",
                      substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("RW", "Reflected Weibull"),
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
                   dldm <- 1/mu - (-y)^sigma 
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   A    <- mu * (-y)^sigma * log(-y)
                   dldd <- 1/sigma + log(-y) - A
                   dldd
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dldm   <- 1/mu - (-y)^sigma 
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   dldm    <- 1/mu - (-y)^sigma
                   A       <- mu * (-y)^sigma * log(-y)
                   dldd    <- 1/sigma + log(-y) - A
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldd2 = function(y, mu, sigma) {
                   A      <- mu * (-y)^sigma * log(-y)
                   dldd   <- 1/sigma + log(-y) - A
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                    G.dev.incr = function(y, mu, sigma, ...) -2*dRW(y, mu, sigma, log=TRUE), 
                         rqres = expression(rqres(pfun="pRW", type="Continuous", y=y, mu=mu, sigma=sigma)), 
                 
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 
                      mu.valid = function(mu) all(mu > 0), 
                   sigma.valid = function(sigma) all(sigma > 0), 
                 
                       y.valid = function(y) all(y < 0)
    ), 
    class=c("gamlss.family", "family"))
}
