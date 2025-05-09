#' The Cosine Sine Exponential family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Cosine Sine Exponential family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dCS2e}
#' 
#' @details 
#' The Cosine Sine Exponential distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\frac{\pi \sigma \mu \exp(\frac{-x} {\nu})}{2 \nu [(\mu\sin(\frac{\pi}{2} \exp(\frac{-x} {\nu})) + \sigma\cos(\frac{\pi}{2} \exp(\frac{-x} {\nu}))]^2}, }
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a CS2e distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_CS2e.R
#'
#' @references
#' Chesneau, C., Bakouch, H. S., & Hussain, T. (2019). A new 
#' class of probability distributions via cosine and sine 
#' functions with applications. Communications in 
#' Statistics-Simulation and Computation, 48(8), 2287-2300.
#' 
#' @importFrom gamlss rqres.plot
#' @export
CS2e <- function (mu.link="log", sigma.link="log", nu.link="log"){
  mstats <- checklink("mu.link", "Cosine Sine Exponential", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Cosine Sine Exponential",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Cosine Sine Exponential", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("CS2e", "Cosine Sine Exponential"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
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
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu))
                   dldm <- 1/mu - ((2*exp1)/ (mu*exp1 + sigma*exp2))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldd <- 1/sigma - ((2*exp2)/ (mu*exp1 + sigma*exp2))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu){
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu))
                   dldv <- y/nu^2 - 1/nu - ((pi*mu* exp(-y/nu)*y*exp2)
                           - (pi*sigma*exp(-y/nu)*y*exp1))/ (nu^2*(mu*exp1 + sigma*exp2))
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu))
                   dldm <- 1/mu - ((2*exp1)/ (mu*exp1 + sigma*exp2))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldm <- 1/mu - ((2*exp1)/ (mu*exp1 + sigma*exp2))
                   dldd <- 1/sigma - ((2*exp2)/ (mu*exp1 + sigma*exp2))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldm <- 1/mu - ((2*exp1)/ (mu*exp1 + sigma*exp2))
                   dldv <- y/nu^2 - 1/nu - ((pi*mu* exp(-y/nu)*y*exp2)
                                            - (pi*sigma*exp(-y/nu)*y*exp1))/ (nu^2*(mu*exp1 + sigma*exp2))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldd <- 1/sigma - ((2*exp2)/ (mu*exp1 + sigma*exp2))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldd <- 1/sigma - ((2*exp2)/ (mu*exp1 + sigma*exp2))
                   dldv <- y/nu^2 - 1/nu - ((pi*mu* exp(-y/nu)*y*exp2)
                                            - (pi*sigma*exp(-y/nu)*y*exp1))/ (nu^2*(mu*exp1 + sigma*exp2))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- sin((pi/2) *exp(-y/nu))
                   exp2 <- cos((pi/2) *exp(-y/nu)) 
                   dldv <- y/nu^2 - 1/nu - ((pi*mu* exp(-y/nu)*y*exp2)
                                            - (pi*sigma*exp(-y/nu)*y*exp1))/ (nu^2*(mu*exp1 + sigma*exp2))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dCS2e(y, mu, sigma, nu, log=TRUE), 
                 rqres      = expression(rqres(pfun="pCS2e", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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
