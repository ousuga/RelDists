#' The Weigted Generalized Exponential-Exponential family
#' 
#' @author Johan David Marin Benjumea, \email{johand.marin@@udea.edu.co}
#' 
#' @description 
#' The Weigted Generalized Exponential-Exponential family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dWGEE}
#' 
#' @details 
#' The Weigted Generalized Exponential-Exponential distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \sigma \nu \exp(-\nu x) (1 - \exp(-\nu x))^{\sigma - 1} (1 - \exp(-\mu \nu x)) / 1 - \sigma B(\mu + 1, \sigma),}
#' 
#' for \eqn{x > 0}, \eqn{\mu > 0}, \eqn{\sigma > 0} and \eqn{\nu > 0}.  
#' 
#' @returns Returns a gamlss.family object which can be used to fit a WGEE distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_WGEE.R 
#' 
#' @references
#'\insertRef{mahdavi2015two}{RelDists}
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
WGEE <- function (mu.link="log", sigma.link="log", nu.link="log"){
  mstats <- checklink("mu.link", "Weigted Generalized Exponential-Exponential", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Weigted Generalized Exponential-Exponential",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Weigted Generalized Exponential-Exponential", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("WGEE", "Weigted Generalized Exponential-Exponential"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
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
                 
                 dldm = function(y, mu, sigma, nu){
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1,sigma)
                   exp3 <- (digamma(mu+1) - digamma(mu+sigma+1))
                   dldm <- (nu*y*exp1^mu)/(1 - exp1^mu) + 
                     (sigma*exp2*exp3)/(1 - sigma*exp2)
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1, sigma)
                   exp4 <- (digamma(sigma) - digamma(mu+sigma+1))
                   dldd <- (1 + sigma*exp4)*exp2/(1 - sigma*exp2) + 1/sigma +
                     log(1 - exp1)  
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu){
                   exp1 <- exp(-nu*y)
                   dldv <- 1/nu - y + (sigma - 1)*((y*exp1)/(1 - exp1)) + 
                     ((mu*y*exp1^mu)/(1 - exp1^mu))
                   dldv
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1,sigma)
                   exp3 <- (digamma(mu+1) - digamma(mu+sigma+1))
                   dldm <- (nu*y*exp1^mu)/(1 - exp1^mu) + 
                     (sigma*exp2*exp3)/(1 - sigma*exp2)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1,sigma)
                   exp3 <- (digamma(mu+1) - digamma(mu+sigma+1))
                   exp4 <- (digamma(sigma) - digamma(mu+sigma+1))
                   dldm <- (nu*y*exp1^mu)/(1 - exp1^mu) + 
                     (sigma*exp2*exp3)/(1 - sigma*exp2)
                   dldd <- (1 + sigma*exp4)*exp2/(1 - sigma*exp2) + 1/sigma +
                     log(1 - exp1)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1,sigma)
                   exp3 <- (digamma(mu+1) - digamma(mu+sigma+1))
                   dldm <- (nu*y*exp1^mu)/(1 - exp1^mu) + 
                     (sigma*exp2*exp3)/(1 - sigma*exp2)
                   dldv <- 1/nu - y + (sigma - 1)*((y*exp1)/(1 - exp1)) + 
                     ((mu*y*exp1^mu)/(1 - exp1^mu))
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1, sigma)
                   exp4 <- (digamma(sigma) - digamma(mu+sigma+1))
                   dldd <- (1 + sigma*exp4)*exp2/(1 - sigma*exp2) + 1/sigma +
                     log(1 - exp1)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   exp2 <- beta(mu+1, sigma)
                   exp4 <- (digamma(sigma) - digamma(mu+sigma+1))
                   dldd <- (1 + sigma*exp4)*exp2/(1 - sigma*exp2) + 1/sigma +
                     log(1 - exp1)
                   dldv <- 1/nu - y + (sigma - 1)*((y*exp1)/(1 - exp1)) + 
                     ((mu*y*exp1^mu)/(1 - exp1^mu))
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   exp1 <- exp(-nu*y)
                   dldv <- 1/nu - y + (sigma - 1)*((y*exp1)/(1 - exp1)) + 
                     ((mu*y*exp1^mu)/(1 - exp1^mu))
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dWGEE(y, mu, sigma, nu,log=TRUE), 
                 rqres = expression(rqres(pfun="pWGEE", type="Continuous", 
                                          y=y, mu=mu, sigma=sigma, nu=nu)),
                 
                 mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu    <- rep(1, length(y))),
                 
                 mu.valid = function(mu) all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu) all(nu > 0), 
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}