#' The Wald family
#'
#' @author Sofia Cuartas García, \email{scuartasg@unal.edu.co}
#'
#' @description
#' The function \code{WALD()} defines the wALD distribution, two-parameter
#' continuous distribution for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma parameter.
#'
#' @references
#' Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to 
#' response time data: An example using functions for the S-PLUS package. 
#' Behavior Research Methods, Instruments, & Computers, 36, 678-694.
#'
#' @seealso \link{dWALD}.
#'
#' @details
#' The Wald distribution with parameters \eqn{\mu} and \eqn{sigma} has density given by
#'
#' \eqn{\operatorname{f}(x |\mu, \sigma)=\frac{\sigma}{\sqrt{2 \pi x^3}} \exp \left[-\frac{(\sigma-\mu x)^2}{2x}\right ], x>0 }
#'
#' @returns Returns a gamlss.family object which can be used to fit a WALD distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_WALD.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
WALD <- function (mu.link="log", sigma.link="log"){
  mstats <- checklink("mu.link", "WALD", 
                      substitute(mu.link), 
                      c("log", "own"))
  dstats <- checklink("sigma.link", "WALD", 
                      substitute(sigma.link), 
                      c("log", "own"))
  
  structure(list(family=c("WALD", "Wald"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Continuous",
                 
                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 
                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 
                 # First derivates
                 dldm = function(y, mu, sigma) {
                   dldm <- sigma-mu*y
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dldd <- 1/sigma - (sigma-mu*y)/y
                   dldd
                 },
                 
                 # Second derivates
                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dWALD(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.0001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dWALD(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.0001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dWALD(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.0001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dWALD(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.0001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dWALD(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pWALD", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial    = expression(mu    <- rep(wald_start(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(wald_start(y)[2], length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 
                 y.valid = function(y) all(y > 0),
                 
                 mean = function(mu, sigma) sigma/mu,
                 variance = function(mu, sigma) sigma/(mu^3)
  ),
  class=c("gamlss.family", "family"))
}
#' Initial values for WALD
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with starting values.
#' @keywords internal
#' @export
#' @importFrom stats var
wald_start <- function(y) {
  mu <- sqrt(mean(y)/var(y))
  sigma <- mu*mean(y)
  return(c(mu, sigma))
}

