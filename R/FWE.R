#' The Flexible Weibull Extension family
#' 
#' @description 
#' The function \code{FWE()} defines the Flexible Weibull distribution, a two parameter 
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @details 
#' The Flexible Weibull extension with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x) = (\mu + \sigma/x^2) exp(\mu x - \sigma/x) exp(-exp(\mu x-\sigma/x))}
#' 
#' for x>0.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a FWE distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_FWE.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
FWE <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "Flexible Weibull Extension", substitute(mu.link), c("log", "identity"))
  dstats <- checklink("sigma.link", "Flexible Weibull Extension", substitute(sigma.link), c("log", "identity"))
  
  structure(list(family = c("FWE", "Flexible Weibull Extension"),
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
                   dldm <- (1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y))
                   dldm
                  },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dldm   <- 1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y)
                   d2ldm2 <- -dldm * dldm
                 },
                 
                 dldd = function(y, mu, sigma){
                   dldd <- (1/(mu*y^2+sigma)-1/y+exp(mu*y-sigma/y)/y) 
                   dldd
                 } ,
                 
                 d2ldd2 = function(y,mu,sigma) {
                   dldd   <- 1/(mu*y^2+sigma) - 1/y + exp(mu*y-sigma/y)/y
                   d2ldd2 <- -dldd * dldd
                 },
                 
                 d2ldmdd = function(y,mu,sigma) {
                   dldm    <- 1/(mu+sigma/y^2) + y - y*exp(mu*y-sigma/y)
                   dldd    <- 1/(mu*y^2+sigma) - 1/y + exp(mu*y-sigma/y)/y
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                    G.dev.incr = function(y, mu, sigma, ...) -2*dFWE(y, mu, sigma, log=TRUE), 
                         rqres = expression(rqres(pfun="pFWE", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial    = expression(    mu <- rep(initValuesFWE(y)[1], length(y)) ),     
                 sigma.initial = expression( sigma <- rep(initValuesFWE(y)[2], length(y)) ), 
                 
                      mu.valid = function(mu) all(mu > 0) , 
                   sigma.valid = function(sigma)  all(sigma > 0), 
                 
                       y.valid = function(y)  all(y > 0)
  ),
  class = c("gamlss.family","family"))
}
#'
#' initValuesFWE
#' 
#' This function generates initial values for FWE distribution.
#' 
#' @param y vector with the random sample
#' @keywords internal
#' @export
#' @importFrom stats coef ecdf lm
initValuesFWE <- function(y) {
  F_hat <- ecdf(y)
  p <- F_hat(y)
  p[p == 1] <- 0.999
  yy <- log(-log(1-p))
  x1 <- y
  x2 <- -1/y
  mod <- lm(yy ~ 0 + x1 + x2)
  res <- coef(mod)
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
