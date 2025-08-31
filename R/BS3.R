#' The Birnbaum-Saunders family - Bourguignon & Gallardo (2022)
#' 
#' @description 
#' The function \code{BS3()} defines The Birnbaum-Saunders, 
#' a two parameter distribution, for a \code{gamlss.family} object 
#' to be used in GAMLSS fitting 
#' using the function \code{gamlss()}.
#' 
#' @param mu.link defines the mu.link, with "log" link as the default 
#' for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default 
#' for the sigma.
#' 
#' @references
#' Bourguignon, M., & Gallardo, D. I. (2022). A new look at the 
#' Birnbaum–Saunders regression model. Applied Stochastic Models in 
#' Business and Industry, 38(6), 935-951.
#' 
#' @seealso \link{dBS3}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{(1-\sigma)y+\mu}{2\sqrt{2\pi\mu\sigma(1-\sigma)}y^{3/2}} \exp{\left[ \frac{-1}{2\sigma} \left( \frac{(1-\sigma)y}{\mu} + \frac{\mu}{(1-\sigma)y} - 2 \right) \right]} }
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{0<\sigma<1}. In this 
#' parameterization 
#' \eqn{Mode(X)=\mu} and 
#' \eqn{Var(X)=(\mu\sigma)^2(1+5\sigma^2/4)}.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a 
#' BS3 distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS3.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS3 <- function(mu.link = "log", sigma.link = "logit"){
  mstats <- checklink("mu.link", "BS3", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "BS3", substitute(sigma.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  structure(
    list(family = c("BS3", "Birnbaum-Saunders - third parameterization"),
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
           term1 <- 1 / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * mu)
           term3 <- (1 / (2 * sigma)) * (((1 - sigma) * y) / (mu^2) - 1 / ((1 - sigma) * y))
           result <- term1 + term2 + term3
           return(result)
         },
         
         dldd = function(y, sigma, mu) {
           A <- ((1 - sigma) * y / mu) + (mu / ((1 - sigma) * y)) - 2
           Aprime <- (-y / mu) + (mu / ((1 - sigma)^2 * y))
           term1 <- -y / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * sigma)
           term3 <- 1 / (2 * (1 - sigma))
           term4 <- -(1 / (2 * sigma)) * Aprime
           term5 <- (1 / (2 * sigma^2)) * A
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         # Second derivatives
         d2ldm2 = function(y,sigma,mu) {
           term1 <- 1 / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * mu)
           term3 <- (1 / (2 * sigma)) * (((1 - sigma) * y) / (mu^2) - 1 / ((1 - sigma) * y))
           result <- term1 + term2 + term3
           return(-result*result)
         },
         
         d2ldd2 = function(y, sigma, mu) {
           A <- ((1 - sigma) * y / mu) + (mu / ((1 - sigma) * y)) - 2
           Aprime <- (-y / mu) + (mu / ((1 - sigma)^2 * y))
           term1 <- -y / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * sigma)
           term3 <- 1 / (2 * (1 - sigma))
           term4 <- -(1 / (2 * sigma)) * Aprime
           term5 <- (1 / (2 * sigma^2)) * A
           result <- term1 + term2 + term3 + term4 + term5
           return(-result*result)
         },
         
         d2ldmdd = function(y,sigma,mu){
           
           term1 <- 1 / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * mu)
           term3 <- (1 / (2 * sigma)) * (((1 - sigma) * y) / (mu^2) - 1 / ((1 - sigma) * y))
           dldm <- term1 + term2 + term3
           
           A <- ((1 - sigma) * y / mu) + (mu / ((1 - sigma) * y)) - 2
           Aprime <- (-y / mu) + (mu / ((1 - sigma)^2 * y))
           term1 <- -y / ((1 - sigma) * y + mu)
           term2 <- -1 / (2 * sigma)
           term3 <- 1 / (2 * (1 - sigma))
           term4 <- -(1 / (2 * sigma)) * Aprime
           term5 <- (1 / (2 * sigma^2)) * A
           dldd <- term1 + term2 + term3 + term4 + term5
           
           d2ldmdd = -dldm * dldd
           d2ldmdd
         },
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS3(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS3", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(0.5, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}