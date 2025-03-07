#' The Marshall-Olkin Extended Weibull family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
#' 
#' @description 
#' The Marshall-Olkin Extended Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dMOEW}
#' 
#' @details 
#' The Marshall-Olkin Extended Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu (\nu x)^{\sigma - 1} exp\{{-(\nu x )^{\sigma}}\}}{\{1-(1-\mu) exp\{{-(\nu x )^{\sigma}}\} \}^{2}},}
#' 
#' for x > 0. 
#' 
#' @returns Returns a gamlss.family object which can be used to fit a MOEW distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_MOEW.R 
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{ghitany2005}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
MOEW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Marshall-Olkin Extended Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall-Olkin Extended Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall-Olkin Extended Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("MOEW", "Marshall-Olkin Extended Weibull"), 
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
                 
                 # Primeras derivadas ---------------------------------
                 dldm = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   dldv
                 },
                 
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   A    <- exp(-(nu * y)^sigma)
                   B    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldm <-  1 / mu - 2 * A / B
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   C    <- (1 / sigma) - log(nu * y) * (nu * y)^sigma + log(nu*y)
                   D    <- (1 - mu) * exp(-(nu * y)^sigma) * log(nu * y) * (nu * y)^sigma
                   E    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldd <- C + (-2 * D / E)
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   G    <- 1 / nu + (sigma - 1) / nu - sigma * y * (nu * y)^(sigma - 1)
                   H    <- - 2 * (1 - mu) * sigma * exp(-(nu * y)^sigma) * y * (nu * y)^(sigma - 1)
                   I    <- 1 - (1 - mu) * exp(-(nu * y)^sigma)
                   dldv <- G + H / I
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMOEW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pMOEW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression(mu       <- rep(estim_mu_sigma_MOEW(y)[1], length(y))), 
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_MOEW(y)[2], length(y))), 
                 nu.initial = expression(nu       <- rep(estim_mu_sigma_MOEW(y)[3], length(y))),
                 
                 mu.valid = function(mu)       all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu)       all(nu > 0),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
#' logLik function for MOEW
#' @description Calculates logLik for MOEW distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @export
logLik_MOEW <- function(logparam=c(0, 0, 0), x){
  return(sum(dMOEW(x     = x,
                   mu    = exp(logparam[1]),
                   sigma = exp(logparam[2]),
                   nu    = exp(logparam[3]),
                   log=TRUE)))
}
#' Initial values for MOEW
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_MOEW <- function(y) {
  mod <- optim(par=c(0, 0, 0),
               fn=logLik_MOEW,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]),
           nu_hat    = exp(mod$par[3]))
  names(res) <- c("mu_hat", "sigma_hat", "nu_hat")
  return(res)
}
