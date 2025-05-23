#' The Generalized Gompertz  family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Generalized Gompertz  family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter. 
#' 
#' @seealso \link{dGGD}
#' 
#' @details
#' The Generalized Gompertz  Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \nu \mu \exp(-\frac{\mu}{\sigma}(\exp(\sigma x - 1))) (1 - \exp(-\frac{\mu}{\sigma}(\exp(\sigma x - 1))))^{(\nu - 1)} ,}
#' 
#' for \eqn{x \geq 0}, \eqn{\mu > 0}, \eqn{\sigma \geq 0} and \eqn{\nu > 0}
#' 
#' @returns Returns a gamlss.family object which can be used to fit a GGD distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_GGD.R
#' 
#' @references
#' #' El-Gohary, A., Alshamrani, A., & Al-Otaibi, A. N. (2013). 
#' The generalized Gompertz distribution. Applied mathematical 
#' modelling, 37(1-2), 13-24.
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GGD <- function (mu.link="log", sigma.link="log", nu.link="log"){
  mstats <- checklink("mu.link", "Generalized Gompertz ", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Gompertz ",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Generalized Gompertz ", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("GGD", "Generalized Gompertz "), 
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
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldm <- 1/mu - 1/sigma*(exp1 - 1) + (nu -1)*exp2*(exp1 - 1)/
                     (sigma*(1-exp2))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldd <- y + mu/sigma^2*(exp1 - 1) - mu/sigma*y*exp1 + (nu -1)* mu*exp2*
                     (sigma*y*exp1 - exp1 + 1)/(sigma^2*(1 - exp2))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu){
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldv  <- 1/nu + log(1 - exp2)
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldm  <- 1/mu - 1/sigma*(exp1 - 1) + (nu -1)*exp2*(exp1 - 1)/
                     (sigma*(1-exp2))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldm <- 1/mu - 1/sigma*(exp1 - 1) + (nu -1)*exp2*(exp1 - 1)/
                     (sigma*(1-exp2))
                   dldd <- y + mu/sigma^2*(exp1 - 1) - mu/sigma*y*exp1 + (nu -1)* mu*exp2*
                     (sigma*y*exp1 - exp1 + 1)/(sigma^2*(1 - exp2))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldm <- 1/mu - 1/sigma*(exp1 - 1) + (nu -1)*exp2*(exp1 - 1)/
                     (sigma*(1-exp2))
                   dldv  <- 1/nu + log(1 - exp2)
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldd <- y + mu/sigma^2*(exp1 - 1) - mu/sigma*y*exp1 + (nu -1)* mu*exp2*
                     (sigma*y*exp1 - exp1 + 1)/(sigma^2*(1 - exp2))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldd <- y + mu/sigma^2*(exp1 - 1) - mu/sigma*y*exp1 + (nu -1)* mu*exp2*
                     (sigma*y*exp1 - exp1 + 1)/(sigma^2*(1 - exp2))
                   dldv <- 1/nu + log(1 - exp2)
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(sigma*y)
                   exp2 <- exp(-mu/sigma*(exp1 - 1))
                   dldv <- 1/nu + log(1 - exp2)
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
       
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dGGD(y, mu, sigma, nu, log=TRUE), 
                 rqres      = expression(rqres(pfun="pGGD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial    = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial    = expression(nu    <- rep(1, length(y))),
                 
                 mu.valid    = function(mu)    all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma >= 0), 
                 nu.valid    = function(nu)    all(nu > 0), 
                 
                 y.valid = function(y) all(y >= 0)
  ), 
  class=c("gamlss.family", "family"))
}




