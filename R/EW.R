#' The Exponentiated Weibull family
#' 
#' @description 
#' The Exponentiated Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @details 
#' The Exponentiated Weibull Distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\nu \mu \sigma x^{\sigma-1} \exp(-\mu x^\sigma) (1-\exp(-\mu x^\sigma))^{\nu-1},}
#' 
#' for x > 0. 
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EW <- function (mu.link="log", sigma.link="log", nu.link="log") 
{
  mstats <- checklink("mu.link",    "Exponentiated Weibull", 
                      substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Exponentiated Weibull", 
                      substitute(nu.link),    c("log", "own"))
  
  structure(list(family=c("EW", "Exponentiated Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
                 mu.link   =as.character(substitute(mu.link)), 
                 sigma.link=as.character(substitute(sigma.link)), 
                 nu.link   =as.character(substitute(nu.link)), 
                 
                 mu.linkfun   =mstats$linkfun, 
                 sigma.linkfun=dstats$linkfun, 
                 nu.linkfun   =vstats$linkfun, 
                 
                 mu.linkinv   =mstats$linkinv, 
                 sigma.linkinv=dstats$linkinv, 
                 nu.linkinv   =vstats$linkinv, 
                 
                 mu.dr   =mstats$mu.eta, 
                 sigma.dr=dstats$mu.eta, 
                 nu.dr   =vstats$mu.eta, 
                 
                 
                 dldm=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- (y^sigma)
                   dexp2dm <- exp(-exp1)*dexp1dm
                   dldm <- 1/mu - y^sigma + ((nu-1)*dexp2dm)/exp2
                   dldm
                 },
                 
                 d2ldm2=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- (y^sigma)
                   dexp2dm <-  exp(-exp1)*dexp1dm
                   d2exp2dm2 <- -exp(-exp1)*(dexp1dm)^2
                   d2ldm2 <- -(-(1/mu^2) + ((nu-1)/exp2^2)*(exp2*d2exp2dm2-dexp2dm^2))^2
                   d2ldm2 
                 }, 
                 
                 dldd=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   dldd <- 1/sigma + log(y) -exp1*log(y) + ((nu-1)*dexp2dd)/exp2
                   dldd
                 },
                 
                 d2ldd2=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1) 
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2exp1dd2 <- log(y)*dexp1dd   
                   d2exp2dd2 <- exp(-exp1)*(d2exp1dd2-dexp1dd^2)
                   d2ldd2 <- -(-(1/sigma^2) -log(y)*dexp1dd - ((nu-1)/exp2^2)*(exp2*d2exp2dd2-dexp2dd^2))^2
                   d2ldd2
                 }, 
                 
                 dldv=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dldv <- 1/nu + log(exp2) 
                   dldv 
                 }, 
                 
                 d2ldv2=function(nu) -(-(1/nu^2))^2, 
                 
                 d2ldmdd=function(y, mu, sigma, nu) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2exp1dmdd <- (y^sigma)*log(y)
                   d2exp2dmdd <- exp(-exp1)*(-dexp1dd*dexp1dm + d2exp1dmdd)
                   d2ldmdd <- -((y^sigma)*log(y) + 
                                  ((nu-1)/exp2^2)*(exp2*d2exp2dmdd - dexp2dm*dexp2dd))^2
                   d2ldmdd
                 },  
                 
                 d2ldmdv=function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dm <- y^sigma
                   dexp2dm <- exp(-exp1)*dexp1dm
                   d2ldmdv <- -(dexp2dm/exp2)^2
                   d2ldmdv
                 },   
                 
                 d2ldddv=function(y, mu, sigma) {
                   exp1 <- mu*(y^sigma)
                   exp2 <- 1-exp(-exp1)
                   dexp1dd <- exp1*log(y)
                   dexp2dd <- exp(-exp1)*dexp1dd
                   d2ldddv <- -(dexp2dd/exp2)^2
                   d2ldddv 
                 },   
                 
                 
                 G.dev.incr=function(y, mu, sigma, nu, ...) -2*dEW(y, mu, sigma, nu, log=TRUE), 
                 rqres=expression(rqres(pfun="pEW", type="Continuous",
                                        y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial=expression( mu <-  rep(1, length(y)) ), 
                 sigma.initial=expression( sigma <- rep(1, length(y)) ), 
                 nu.initial=expression( nu <- rep(1, length(y)) ), 
                 
                 mu.valid=function(mu) all(mu >  0), 
                 sigma.valid=function(sigma) all(sigma >  0), 
                 nu.valid=function(nu) all(nu > 0), 
                 
                 y.valid=function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}