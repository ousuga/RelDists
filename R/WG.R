#' The Weibull Geometric family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Weibull Geometric distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dWG}
#' 
#' @details 
#' The weibull geometric distribution with parameters \code{mu},
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = (\sigma \mu^\sigma (1-\nu) x^(\sigma - 1) \exp(-(\mu x)^\sigma)) 
#' (1- \nu \exp(-(\mu x)^\sigma))^{-2},}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{0 < \nu < 1}.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a WG distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_WG.R
#' 
#' @references
#' Barreto-Souza, W., de Morais, A. L., & Cordeiro, G. M. (2011). 
#' The Weibull-geometric distribution. Journal of Statistical 
#' Computation and Simulation, 81(5), 645-657.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
WG <- function (mu.link = "log", sigma.link = "log", nu.link = "logit") {
  mstats <- checklink("mu.link",    "Weibull Geometric", substitute(mu.link),    c("identity", "own"))
  dstats <- checklink("sigma.link", "Weibull Geometric", substitute(sigma.link), c("identity", "own"))
  vstats <- checklink("nu.link",    "Weibull Geometric", substitute(nu.link),    c("logit","probit", "own"))
  
  structure(list(family = c("WG", "Weibull Geometric"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                 nopar = 3, 
                 type = "Continuous", 
                 
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
                   exp1 <- exp(-(mu*y)^sigma)
                   exp2 <- 2/(1 - nu*exp1)
                   dldm <- sigma/mu - sigma*y*(mu*y)^(sigma - 1) -
                     exp2*nu*exp1*sigma*(mu*y)^(sigma-1)*y
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- exp(-(mu*y)^sigma)
                   exp2 <- 2/(1 - nu*exp1)
                   dldd <- 1/sigma + log(mu) + log(y) - (mu*y)^sigma*log(mu*y) -
                     exp2*nu*exp1*(mu*y)^sigma*log(mu*y)
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   exp1 <- exp(-(mu*y)^sigma)
                   exp2 <- 2/(1 - nu*exp1)
                   dldv <- -(1/(1-nu)) + exp2*exp1
                   dldv 
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   exp1   <- exp(-(mu*y)^sigma)
                   exp2   <- 2/(1 - nu*exp1)
                   dldm   <- sigma/mu - sigma*y*(mu*y)^(sigma - 1) -
                     exp2*nu*exp1*sigma*(mu*y)^(sigma-1)*y
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   exp1    <- exp(-(mu*y)^sigma)
                   exp2    <- 2/(1 - nu*exp1)
                   dldm    <- sigma/mu - sigma*y*(mu*y)^(sigma - 1) -
                     exp2*nu*exp1*sigma*(mu*y)^(sigma-1)*y
                   dldd    <- 1/sigma + log(mu) + log(y) - (mu*y)^sigma*log(mu*y) -
                     exp2*nu*exp1*(mu*y)^sigma*log(mu*y)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   exp1    <- exp(-(mu*y)^sigma)
                   exp2    <- 2/(1 - nu*exp1)
                   dldm    <- sigma/mu - sigma*y*(mu*y)^(sigma - 1) -
                     exp2*nu*exp1*sigma*(mu*y)^(sigma-1)*y
                   dldv    <- -(1/(1-nu)) + exp2*exp1
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   exp1   <- exp(-(mu*y)^sigma)
                   exp2   <- 2/(1 - nu*exp1)
                   dldd   <- 1/sigma + log(mu) + log(y) - (mu*y)^sigma*log(mu*y) -
                     exp2*nu*exp1*(mu*y)^sigma*log(mu*y)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   exp1    <- exp(-(mu*y)^sigma)
                   exp2    <- 2/(1 - nu*exp1)
                   dldd    <- 1/sigma + log(mu) + log(y) - (mu*y)^sigma*log(mu*y) -
                     exp2*nu*exp1*(mu*y)^sigma*log(mu*y)
                   dldv    <- -(1/(1-nu)) + exp2*exp1
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 =  function(y, mu, sigma, nu, tau) {
                   exp1   <- exp(-(mu*y)^sigma)
                   exp2   <- 2/(1 - nu*exp1)
                   dldv   <- -(1/(1-nu)) + exp2*exp1
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                    G.dev.incr = function(y, mu, sigma, nu, ...) -2*dWG(y, mu, sigma, nu, log = TRUE), 
                         rqres = expression(rqres(pfun = "pWG", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)), 
                 
                    mu.initial = expression( mu <-  rep(0.5, length(y)) ), 
                 sigma.initial = expression( sigma <- rep(0.5, length(y)) ), 
                    nu.initial = expression( nu <- rep(0.5, length(y)) ), 
                 
                      mu.valid = function(mu) all(mu >  0), 
                   sigma.valid = function(sigma) all(sigma >  0), 
                      nu.valid = function(nu) all(nu > 0 & nu < 1), 
                 
                       y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}