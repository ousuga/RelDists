#' The Gamma Weibull family
#' 
#' @description 
#' The Gamma Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dGammaW}
#' 
#' @details 
#' The Gamma Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)= \frac{\sigma \mu^{\nu}}{\Gamma (\nu)} x^{\nu \sigma - 1} \exp(-\mu x^\sigma),}
#' 
#' for x > 0. 
#' 
#' @examples #' 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rGammaW(n=100, mu = 0.5, sigma = 2, nu=1)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='GammaW',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma and nu
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' exp(coef(mod, what='nu'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n     <- 200
#' x1    <- runif(n)
#' x2    <- runif(n)
#' mu    <- exp(-1.6 * x1)
#' sigma <- exp(1.1 - 1 * x2)
#' nu    <- 1
#' x     <- rGammaW(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, mu.fo=~x1, sigma.fo=~x2, nu.fo=~1, family=GammaW,
#'               control=gamlss.control(n.cyc=50000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' coef(mod, what='nu')
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GammaW <- function (mu.link="log", sigma.link="log", nu.link="log") 
{
  mstats <- checklink("mu.link", "Gamma Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Gamma Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Gamma Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("GammaW", "Gamma Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                 nopar=3, 
                 type="Continuous", 
                 
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
                   dldm <- (nu/mu) - y^sigma
                   dldm
                 },
                 
                 dldd    =  function(y, mu, sigma, nu) {
                   dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                   dldd
                 },
                 
                 dldv    =  function(y, mu, sigma, nu){
                   dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                 dldv
                 },
                 
                d2ldm2   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  d2ldm2 <- -dldm * dldm
                  d2ldm2
                  },
                
                d2ldmdd   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  d2ldmdd <- -dldm * dldd
                  d2ldmdd
                  },
                
                d2ldmdv   =  function(y, mu, sigma, nu, tau) {
                  dldm <- (nu/mu) - y^sigma
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldmdv <- -dldm * dldv
                  d2ldmdv
                  },
                
                d2ldd2   =  function(y, mu, sigma, nu, tau) {
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  d2ldd2 <- -dldd * dldd
                  d2ldd2
                  },
                
                d2ldddv   =  function(y, mu, sigma, nu, tau) {
                  dldd  <- 1/sigma + nu*log(y) - mu*log(y)*y^sigma
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldddv <- -dldd * dldv
                  d2ldddv
                  },
                
                d2ldv2   =  function(y, mu, sigma, nu, tau) {
                  dldv  <- log(mu) + sigma*log(y) - base::digamma(nu)
                  d2ldv2 <- -dldv * dldv
                  d2ldv2
                  },
                
                
                G.dev.incr = function(y, mu, sigma, nu, ...) -2*dGammaW(y, mu, sigma, nu, log=TRUE), 
                rqres      = expression(rqres(pfun="pGammaW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                
                mu.initial    = expression(mu    <- rep(1, length(y))), 
                sigma.initial = expression(sigma <- rep(1, length(y))), 
                nu.initial    = expression(nu    <- rep(1, length(y))),
                
                mu.valid    = function(mu)    all(mu > 0), 
                sigma.valid = function(sigma) all(sigma > 0), 
                nu.valid    = function(nu)    all(nu > 0),
                
                y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}

