#' The two-parameter Chris-Jerry distribution family
#' 
#' @author Manuel Gutierrez Tangarife, \email{mgutierrezta@unal.edu.co}
#' 
#' @description 
#' The function \code{CJ2()} defines The  two-parameter  Chris-Jerry distribution, 
#' a two parameter distribution, for a \code{gamlss.family} object to be used 
#' in GAMLSS fitting using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for 
#' the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default 
#' for the sigma.
#' 
#' @seealso \link{dCJ2}
#' 
#' @details 
#' The two-parameter  Chris-Jerry distribution with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{
#' f(x; \sigma, \mu) = \frac{\mu^2}{\sigma \mu + 2} (\sigma + \mu x^2) e^{-\mu x}; \quad x > 0, \quad \mu > 0, \quad \sigma > 0
#' }
#' 
#' Note: In this implementation we changed the original parameters \eqn{\theta} for \eqn{\mu}
#' and \eqn{\lambda} for \eqn{\sigma} we did it to implement this distribution
#' within gamlss framework.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a CJ2 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_CJ2.R
#' 
#' @references
#' Chinedu, Eberechukwu Q., et al. "New lifetime distribution with applications 
#' to single acceptance sampling plan and scenarios of increasing hazard 
#' rates" Symmetry 15.10 (2023): 188.
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
CJ2 <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "a two-parameter
 Chris-Jerry distribution", substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "a two-parameter
 Chris-Jerry distribution", substitute(sigma.link), c("log", "identity"))
  
  structure(list(family = c("CJ2", "a two-parameter
 Chris-Jerry distribution"),
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
                   dm   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=1e-8)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                   
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=1e-8)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=1e-8)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },
                 
                 d2ldd2 = function(y,mu,sigma) {
                   dd   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=1e-8)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2 
                 },
                 
                 d2ldmdd = function(y,mu,sigma) {
                   dm   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=1e-8)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dCJ2(y, mu, sigma,log=TRUE),
                                                 theta="sigma",
                                                 delta=1e-8)
                   dldd <- as.vector(attr(dd, "gradient"))
                   
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },
                 
                 G.dev.incr = function(y, mu, sigma, ...) -2*dCJ2(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pCJ2", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial    = expression(    mu <- rep(1, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 y.valid = function(y)  all(y >= 0)
  ),
  class = c("gamlss.family","family"))
}
#' logLik function for CJ2
#' @description Calculates logLik for CJ2 distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @export
logLik_CJ2 <- function(logparam=c(0, 0), x){
  return(sum(dCJ2(x     = x,
                  mu    = exp(logparam[1]),
                  sigma = exp(logparam[2]),
                  log=TRUE)))
}
#'
#' initValuesCJ2
#' 
#' This function generates initial values for CJ2 distribution.
#' 
#' @param y vector with the random sample
#' @keywords internal
#' 
#' @return 
#' A two-length numeric vector with initial estimates for \eqn{mu} and \eqn{sigma} 
#' parameters from CJ2 distribution (see \code{\link{dCJ2}}).
#' 
#' @export
#' @importFrom stats optim
estim_mu_sigma_CJ2 <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_CJ2,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
