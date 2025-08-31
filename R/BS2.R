#' The Birnbaum-Saunders family - Santos-Neto et al. (2014)
#' 
#' @description 
#' The function \code{BS2()} defines The Birnbaum-Saunders, 
#' a two parameter distribution, for a \code{gamlss.family} object 
#' to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default 
#' for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default 
#' for the sigma.
#' 
#' @references
#' Santos-Neto, M., Cysneiros, F. J. A., Leiva, V., & Barros, M. 
#' (2014). A reparameterized Birnbaum–Saunders distribution 
#' and its moments, estimation and applications. 
#' REVSTAT-Statistical Journal, 12(3), 247-272.
#' 
#' @seealso \link{dBS2}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{\exp(\sigma/2)\sqrt{\sigma+1}}{4\sqrt{\pi\mu}x^{3/2}} \left[ x + \frac{\mu\sigma}{\sigma+1} \right] \exp\left( \frac{-\sigma}{4} \left(\frac{x(\sigma+1)}{\mu\sigma}+\frac{\mu\sigma}{x(\sigma+1)} \right) \right) }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization 
#' \eqn{E(X)=\mu} and 
#' \eqn{Var(X)=\frac{\mu^2(2\sigma+5)}{(\sigma+1)^2}}.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a 
#' BS2 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS2.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS2 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS2", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "BS2", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("BS2", "Birnbaum-Saunders - second parameterization"),
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
         
         # First derivatives
         dldm = function(y, sigma, mu) {
           term1 <- -1 / (2 * mu)
           term2 <- (sigma / (sigma + 1)) / (y + (sigma * mu) / (sigma + 1))
           term3 <- y * (sigma + 1) / (4 * mu^2)
           term4 <- -sigma^2 / (4 * y * (sigma + 1))
           result <- term1 + term2 + term3 + term4
           return(result)
         },
         
         dldd = function(y, sigma, mu) {
           term1 <- 1/2
           term2 <- 1 / (2 * (sigma + 1))
           term3 <- mu / ((sigma + 1)^2 * (y + (sigma * mu) / (sigma + 1)))
           term4 <- -y / (4 * mu)
           term5 <- -mu * (sigma^2 + 2 * sigma) / (4 * y * (sigma + 1)^2)
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         # Second derivatives
         d2ldm2 = function(y,sigma,mu) {
           term1 <- -1 / (2 * mu)
           term2 <- (sigma / (sigma + 1)) / (y + (sigma * mu) / (sigma + 1))
           term3 <- y * (sigma + 1) / (4 * mu^2)
           term4 <- -sigma^2 / (4 * y * (sigma + 1))
           result <- term1 + term2 + term3 + term4
           return(-result*result)
         },
         
         d2ldd2 = function(y, sigma, mu) {
           term1 <- 1/2
           term2 <- 1 / (2 * (sigma + 1))
           term3 <- mu / ((sigma + 1)^2 * (y + (sigma * mu) / (sigma + 1)))
           term4 <- -y / (4 * mu)
           term5 <- -mu * (sigma^2 + 2 * sigma) / (4 * y * (sigma + 1)^2)
           result <- term1 + term2 + term3 + term4 + term5
           return(-result*result)
         },
         
         d2ldmdd = function(y,sigma,mu){
           
           term1 <- -1 / (2 * mu)
           term2 <- (sigma / (sigma + 1)) / (y + (sigma * mu) / (sigma + 1))
           term3 <- y * (sigma + 1) / (4 * mu^2)
           term4 <- -sigma^2 / (4 * y * (sigma + 1))
           dldm <- term1 + term2 + term3 + term4
           
           term1 <- 1/2
           term2 <- 1 / (2 * (sigma + 1))
           term3 <- mu / ((sigma + 1)^2 * (y + (sigma * mu) / (sigma + 1)))
           term4 <- -y / (4 * mu)
           term5 <- -mu * (sigma^2 + 2 * sigma) / (4 * y * (sigma + 1)^2)
           dldd <- term1 + term2 + term3 + term4 + term5
           
           d2ldmdd = -dldm * dldd
           d2ldmdd
         },
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS2(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS2", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(1, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}