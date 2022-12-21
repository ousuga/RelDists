#' The Odd Weibull family
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' 
#' @description 
#' The function \code{OW()} defines the Odd Weibull distribution, a three parameter 
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu.
#' 
#' @details 
#' The odd Weibull with parameters \code{mu}, \code{sigma} and \code{nu}
#' has density given by
#' 
#' \eqn{f(t) = \left( \frac{\sigma\nu}{t} \right) (\mu t)^\sigma
#'      e^{(\mu t)^\sigma} \left(e^{(\mu t)^{\sigma}}-1\right)^{\nu-1}
#'      \left[ 1 + \left(e^{(\mu t)^{\sigma}}-1\right)^\nu \right]^{-2}}
#'
#' for x > 0.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a OW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_OW.R
#' 
#' @references
#' \insertRef{Cooray2006}{RelDists}
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
OW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link",    "Odd Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Odd Weibull", substitute(sigma.link), c("identity", "log", "own"))
  vstats <- checklink("nu.link",    "Odd Weibull", substitute(nu.link),    c("identity", "log", "own"))
  
  # valid_values <- OW_modifications(valid.values)
  # sigma.space <- valid_values$sigma.space
  # nu.space <- valid_values$nu.space
  
  structure(list(family = c("OW", "Odd Weibull"),
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
                 nopar = 3,
                 type = "Continuous",
                 
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
                   prod1 <- (mu*y)^sigma
                   expT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1m <- ( sigma*prod1 )/mu
                   T2m <- sigma*prod1*(nu - 1)*exp(prod1)/
                     ( mu*expm1(prod1) )
                   T3m <- -( 2*nu*sigma*prod1*exp(prod1)*(expm1(prod1))^(nu-1) )/
                     (mu*(expT1))
                   dldm <- sigma/mu + T1m + T2m + T3m
                   as.numeric(dldm)
                 },
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expmT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1m <- ( sigma*prod1 )/mu
                   T2m <- sigma*prod1*(nu - 1)*exp(prod1)/
                     ( mu*expm1(prod1) )
                   T3m <- -( 2*nu*sigma*prod1*exp(prod1)*(expm1(prod1))^(nu-1) )/
                     (mu*(expmT1))
                   dldm <- sigma/mu + T1m + T2m + T3m
                   d2ldm2 <- -dldm*dldm
                   as.numeric(d2ldm2)
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expdT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1d <- (nu - 1)*exp(prod1)/expm1(prod1)
                   T2d <- -2*nu*exp(prod1)*(expm1(prod1))^(nu - 1)/expdT1
                   dldd <- 1/sigma + log(mu*y) + log(mu*y)*prod1*(1 + T1d + T2d)
                   as.numeric(dldd)
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expdT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1d <- (nu - 1)*exp(prod1)/expm1(prod1)
                   T2d <- -2*nu*exp(prod1)*(expm1(prod1))^(nu - 1)/expdT1
                   dldd <- 1/sigma + log(mu*y) + log(mu*y)*prod1*(1 + T1d + T2d)
                   d2ldd2 <- -dldd*dldd
                   as.numeric(d2ldd2)
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expvT1 <- 1 + ( expm1( prod1 ) )^nu
                   dldv <- 1/nu + log(expm1(prod1))*(1 - 2*(expm1(prod1))^nu/expvT1)
                   as.numeric(dldv)
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expvT1 <- 1 + ( expm1( prod1 ) )^nu
                   dldv <- 1/nu + log(expm1(prod1))*(1 - 2*(expm1(prod1))^nu/expvT1)
                   d2ldv2 <- -dldv*dldv
                   as.numeric(d2ldv2)
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1m <- ( sigma*prod1 )/mu
                   T2m <- sigma*prod1*(nu - 1)*exp(prod1)/
                     ( mu*expm1(prod1) )
                   T3m <- -( 2*nu*sigma*prod1*exp(prod1)*(expm1(prod1))^(nu-1) )/
                     (mu*(expT1))
                   dldm <- sigma/mu + T1m + T2m + T3m
                   
                   expdT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1d <- (nu - 1)*exp(prod1)/expm1(prod1)
                   T2d <- -2*nu*exp(prod1)*(expm1(prod1))^(nu - 1)/expdT1
                   dldd <- 1/sigma + log(mu*y) + log(mu*y)*prod1*(1 + T1d + T2d)
                   
                   d2ldmdd <- -dldm*dldd
                   as.numeric(d2ldmdd)
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1m <- ( sigma*prod1 )/mu
                   T2m <- sigma*prod1*(nu - 1)*exp(prod1)/
                     ( mu*expm1(prod1) )
                   T3m <- -( 2*nu*sigma*prod1*exp(prod1)*(expm1(prod1))^(nu-1) )/
                     (mu*(expT1))
                   dldm <- sigma/mu + T1m + T2m + T3m
                   
                   expvT1 <- 1 + ( expm1( prod1 ) )^nu
                   dldv <- 1/nu + log(expm1(prod1))*(1 - 2*(expm1(prod1))^nu/expvT1)
                   
                   d2ldmdv <- -dldm*dldv
                   as.numeric(d2ldmdv)
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   prod1 <- (mu*y)^sigma
                   expdT1 <- 1 + ( expm1( prod1 ) )^nu
                   T1d <- (nu - 1)*exp(prod1)/expm1(prod1)
                   T2d <- -2*nu*exp(prod1)*(expm1(prod1))^(nu - 1)/expdT1
                   dldd <- 1/sigma + log(mu*y) + log(mu*y)*prod1*(1 + T1d + T2d)
                   
                   expvT1 <- 1 + ( expm1( prod1 ) )^nu
                   dldv <- 1/nu + log(expm1(prod1))*(1 - 2*(expm1(prod1))^nu/expvT1)
                   
                   d2ldddv <- -dldd*dldv
                   as.numeric(d2ldddv)
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dOW(y, mu, sigma, nu, log = TRUE),
                 rqres = expression(rqres(pfun = "pOW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu)),
                 
                 mu.initial = expression(mu <- rep(1/mean(y), length(y))),
                 
                 # Increasing hazard as default
                 sigma.initial = expression(sigma <- rep(2, length(y))),
                 nu.initial = expression(nu <- rep(6, length(y))),
                 
                 mu.valid = function(mu) all(mu >  0),
                 sigma.valid = function(sigma) all(sigma > 1),
                 nu.valid = function(nu) all(nu >  1),
                 # sigma.valid = sigma.space,
                 # nu.valid = nu.space,
                 
                 y.valid = function(y) all(y > 0)
  ),
  class = c("gamlss.family", "family"))
}
