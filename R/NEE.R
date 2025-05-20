#' New Exponentiated Exponential family
#'
#' @description
#' The function \code{NEE()} defines the New Exponentiated Exponential distribution, a two parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default for the sigma.
#' 
#' @seealso \link{dNEE}
#' 
#' @details
#' The New Exponentiated Exponential distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x | \mu, \sigma) = \log(2^\sigma) \mu \exp(-\mu x) (1-\exp(-\mu x))^{\sigma-1} 2^{(1-\exp(-\mu x))^\sigma}, }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}.
#' 
#' Note: In this implementation we changed the original parameters 
#' \eqn{\theta} for \eqn{\mu} and \eqn{\alpha} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a 
#' NEE distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_NEE.R
#' 
#' @references
#' Hassan, Anwar, I. H. Dar, and M. A. Lone. "A New Class of Probability 
#' Distributions With An Application to Engineering Data." 
#' Pakistan Journal of Statistics and Operation Research 20.2 (2024): 217-231.
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats dnorm pnorm
#' @export
NEE <- function(mu.link="log", sigma.link="log") {
  
  mstats <- checklink("mu.link", "New Exponentiated Exponential",
                      substitute(mu.link), c("log", "inverse", "own"))
  dstats <- checklink("sigma.link", "New Exponentiated Exponential",
                      substitute(sigma.link), c("log", "inverse", "own"))
  
  structure(list(family=c("NEE", "New Exponentiated Exponential"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",
                 
                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
      
                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 
                 # Primeras derivadas
                 
                 dldm = function(y, mu, sigma) {
                   w <- log(2^sigma)
                   z <- exp(-mu*y)
                   r <- 1 - z
                   
                   p1 <- 1/mu 
                   p2 <- y 
                   p3 <- (sigma -1)*y*z / r 
                   p4 <- sigma*y*z*(r^(sigma-1))*log(2) 
                   dldm <- p1 - p2 + p3 + p4
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   w <- log(2^sigma)
                   z <- exp(-mu*y)
                   r <- 1 - z
                   
                   p1 <- log(2)/w
                   p2 <- log(r)
                   p3 <- (r^sigma)*log(r)*log(2)
                   dldd <- p1 + p2 + p3 
                   dldd
                 },
                 
                 # Segundas derivadas
                 
                 d2ldm2 = function(y, mu, sigma) {
                   w <- log(2^sigma)
                   z <- exp(-mu*y)
                   r <- 1 - z
                   
                   p1 <- 1/mu 
                   p2 <- y 
                   p3 <- (sigma -1)*y*z / r 
                   p4 <- sigma*y*z*(r^(sigma-1))*log(2)
                   dldm <- p1 - p2 + p3 + p4
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   w <- log(2^sigma)
                   z <- exp(-mu*y)
                   r <- 1 - z
                   
                   p1 <- log(2)/w
                   p2 <- log(r)
                   p3 <- (r^sigma)*log(r)*log(2)
                   dldd <- p1 + p2 + p3
                   d2ldmdd <- - dldd * dldd
                   d2ldmdd
                 },
                 
                 d2ldd2  = function(y, mu, sigma) {
                   w <- log(2^sigma)
                   z <- exp(-mu*y)
                   r <- 1 - z
                   
                   p1 <- log(2)/w
                   p2 <- log(r)
                   p3 <- (r^sigma)*log(r)*log(2)
                   dldd <- p1 + p2 + p3
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, ...) -2*dNEE(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pNEE", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_NEE(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_NEE(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 
                 y.valid = function(y)  all(y >= 0)
  ),
  class=c("gamlss.family", "family"))
}
#'
#' estim_mu_sigma_NEE
#'
#' This function generates initial values for NEE distribution.
#'
#' @param y vector with the random sample
#' @examples
#' y <- rNEE(n=100, mu=4.6, sigma=0.35)
#' estim_mu_sigma_NEE(y=y)
#' @importFrom stats optim
#' @export
estim_mu_sigma_NEE <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_NEE,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#'
#' logLik_NEE
#'
#' This is an auxiliar function to obtain the logLik for NEE.
#'
#' @param param vector with the values for mu and sigma
#' @param x vector with the data
#' @examples
#' y <- rNEE(n=100, mu=1, sigma=1)
#' logLik_NEE(param=c(0, 0), x=y)
#' @importFrom stats optim
#' @export
logLik_NEE <- function(param=c(0, 0), x){
  return(sum(dNEE(x,
                  mu    = exp(param[1]),
                  sigma = exp(param[2]),
                  log=TRUE)))
}
