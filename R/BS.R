#' The Birnbaum-Saunders family
#' 
#' @description 
#' The function \code{BS()} defines The Birnbaum-Saunders, 
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
#' Birnbaum, Z.W. and Saunders, S.C. (1969a). A new family of life 
#' distributions. J. Appl. Prob., 6, 319-327.
#' 
#' Roquim, F. V., Ramires, T. G., Nakamura, L. R., Righetto, A. J., 
#' Lima, R. R., & Gomes, R. A. (2021). Building flexible regression 
#' models: including the Birnbaum-Saunders distribution in the 
#' gamlss package. Semina: Ciências Exatas e Tecnológicas, 
#' 42(2), 163-168.
#' 
#' @seealso \link{dBS}.
#' 
#' @details 
#' The Birnbaum-Saunders with parameters \code{mu} and \code{sigma}
#' has density given by
#' 
#' \eqn{f(x|\mu,\sigma) = \frac{x^{-3/2}(x+\mu)}{2\sigma\sqrt{2\pi\mu}} \exp\left(\frac{-1}{2\sigma^2}(\frac{x}{\mu}+\frac{\mu}{x}-2)\right)}
#' 
#' for \eqn{x>0}, \eqn{\mu>0} and \eqn{\sigma>0}. In this 
#' parameterization \eqn{\mu} is the median of \eqn{X}, 
#' \eqn{E(X)=\mu(1+\sigma^2/2)} and 
#' \eqn{Var(X)=(\mu\sigma)^2(1+5\sigma^2/4)}. The functions 
#' proposed here
#' corresponds to the functions created by Roquim et al. (2021)
#' with minor modifications to obtain correct log-likelihoods
#' and random samples.
#' 
#' @returns Returns a gamlss.family object which can be used to fit a 
#' BS distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_BS.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
BS <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "BS", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("BS", "Birnbaum-Saunders"),
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
         dldm = function(y,sigma,mu) 1/(y+mu)-1/(2*mu)-1/(2*sigma^2)*(1/y-y/mu^2),
         dldd = function(y,sigma,mu) 1/sigma^3*(y/mu+mu/y-2)-1/sigma,
         
         # Second derivatives
         d2ldm2 = function(y,sigma,mu) -(1/(y+mu)-1/(2*mu)-1/(2*sigma^2)*(1/y-y/mu^2))^2,
         d2ldd2 = function(y,sigma,mu) -2/sigma^2,
         
         d2ldmdd = function(y,sigma,mu){
           nd = gamlss::numeric.deriv(dBS(y, mu, sigma,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd, "gradient"))
           nd = gamlss::numeric.deriv(dBS(y, mu, sigma, log = TRUE), "sigma", delta = 1e-04)
           dldd = as.vector(attr(nd, "gradient"))
           d2ldmdd = -dldm * dldd
           d2ldmdd
         },
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(median(y), length(y))}),
         sigma.initial = expression({sigma <- rep(1, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}