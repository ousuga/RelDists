#' The Modified Weibull family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
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
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu \geq 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a MW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_MW.R 
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Lai, C. D., Xie, M., & Murthy, D. N. P. (2003). 
#' A modified Weibull distribution. 
#' IEEE Transactions on reliability, 52(1), 33-37.
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
