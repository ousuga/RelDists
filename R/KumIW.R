#' The Kumaraswamy Inverse Weibull family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Kumaraswamy Inverse Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dKumIW}
#' 
#' @details 
#' The Kumaraswamy Inverse Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \mu \sigma \nu x^{-\mu - 1} \exp{- \sigma x^{-\mu}} (1 - \exp{- \sigma x^{-\mu}})^{\nu - 1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a KumIW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_KumIW.R
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#'
#' Shahbaz, M. Q., Shahbaz, S., & Butt, N. S. (2012). 
#' The Kumaraswamy Inverse Weibull Distribution. 
#' Pakistan journal of statistics and operation research, 479-489.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
KumIW <- function (mu.link="log", sigma.link="log", nu.link="log"){
  mstats <- checklink("mu.link", "Kumaraswamy Inverse-Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Kumaraswamy Inverse-Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Kumaraswamy Inverse-Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("KumIW", "Kumaraswamy Inverse-Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=4, 
                 type="Continuous", 
                 
                 mu.link    = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link    = as.character(substitute(nu.link)), 
                 
                 mu.linkfun    = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun    = vstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv    = vstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr    = vstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldm <- 1/mu - log(y) + exp1*log(y) - exp3*exp1*log(y)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldd <- 1/sigma - y^(-mu) + exp3*y^(-mu) 
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu){
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   dldv <- 1/nu + log(exp2)
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldm  <-  1/mu - log(y) + exp1*log(y) - exp3*exp1*log(y)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldm <- 1/mu - log(y) + exp1*log(y) - exp3*exp1*log(y)
                   dldd <- 1/sigma - y^(-mu) + exp3*y^(-mu)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldm <-  1/mu - log(y) + exp1*log(y) - exp3*exp1*log(y)
                   dldv <- 1/nu + log(exp2)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldd <- 1/sigma - y^(-mu) + exp3*y^(-mu)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   exp3 <- (nu - 1) * (exp(-exp1)/(exp2))
                   dldd <- 1/sigma - y^(-mu) + exp3*y^(-mu)
                   dldv  <- 1/nu + log(exp2)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- sigma * y^(-mu)
                   exp2 <- 1- exp(-exp1)
                   dldv  <- 1/nu + log(exp2)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dKumIW(y, mu, sigma, nu, log=TRUE), 
                 rqres      = expression(rqres(pfun="pKumIW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial    = expression(nu    <- rep(1, length(y))),
                 
                 mu.valid    = function(mu)    all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid    = function(nu)    all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}




