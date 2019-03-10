#' The Generalized modified Weibull Distribution 
#' 
#' @description 
#' The Generalized modified Weibull Distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "sqrt" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "sqrt" link as the default for the tau parameter. 
#' 
#' @seealso \link{dGMW}
#' 
#' @details 
#' The Generalized modified Weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#'\eqn{f(x)= \mu \sigma x^{\nu - 1}(\nu + \tau x) \exp(\tau x - \mu x^{\nu} e^{\tau x})
#' [1 - \exp(- \mu x^{\nu} e^{\tau x})]^{\sigma-1},}
#' 
#' for x > 0. 
#' 
#' @examples  
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rGMW(n=100, mu=2, sigma=0.5, nu=2, tau=1.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~ 1, family='GMW',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma and nu
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' (coef(mod, what='nu'))^2
#' (coef(mod, what='tau'))^2
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' mu <- exp(2 + -3 * x1)
#' sigma <- exp(3 - 2 * x2)
#' nu <- 2
#' tau <- 1.5
#' x <- rGMW(n=n, mu, sigma, nu, tau)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~ 1, family="GMW", 
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' (coef(mod, what="nu"))^2
#' (coef(mod, what="tau"))^2
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GMW <- function (mu.link = "log", sigma.link = "log", nu.link = "sqrt", tau.link = "sqrt") 
{
  mstats <- checklink("mu.link",    "Generalized Modified Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Modified Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Generalized Modified Weibull", substitute(nu.link),    c("sqrt", "own"))
  tstats <- checklink("tau.link",   "Generalized Modified Weibull", substitute(tau.link),    c("sqrt", "own"))
  
  structure(
   list(family = c("GMW", "Generalized Modified Weibull"), 
    parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE), 
         nopar = 4, 
          type = "Continuous", 
                 
    mu.link    = as.character(substitute(mu.link)), 
    sigma.link = as.character(substitute(sigma.link)), 
    nu.link    = as.character(substitute(nu.link)), 
    tau.link   = as.character(substitute(tau.link)), 
                 
 mu.linkfun    = mstats$linkfun, 
 sigma.linkfun = dstats$linkfun, 
 nu.linkfun    = vstats$linkfun, 
 tau.linkfun   = tstats$linkfun, 
                 
 mu.linkinv    = mstats$linkinv, 
 sigma.linkinv = dstats$linkinv, 
 nu.linkinv    = vstats$linkinv, 
 tau.linkinv   = tstats$linkinv,
                 
      mu.dr    = mstats$mu.eta, 
      sigma.dr = dstats$mu.eta, 
      nu.dr    = vstats$mu.eta, 
      tau.dr   = tstats$mu.eta,

# Inicio de las derivadas

# Primeras derivadas
dldm = function(y, mu, sigma, nu, tau) {
  dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/(1-exp(-mu*(y^nu)*exp(tau*y))))
  dldm
  },
                
dldd = function(y, mu, sigma, nu, tau) {
  dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
  dldd
  },
            
dldv = function(y, mu, sigma, nu, tau) {
  part1 <- mu*exp(tau*y)*(y^nu)*log(y)
  part2 <- mu*(y^nu)*exp(tau*y)
  dldv  <- log(y) + 1/(nu + tau*y) - part1 + (((sigma-1)*part1*exp(-part2))/(1-exp(-part2)))
  dldv
  },
            
dldt = function(y, mu, sigma, nu, tau) {
  part1 <- mu*(y^nu)*exp(tau*y)
  dldt  <- (y/(nu+(tau*y))) + y -y*part1 +(((sigma-1)*exp(-part1)*y*part1)/(1-exp(-part1)))
  dldt
  },
            
# Segundas derivadas
            
d2ldm2 = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "mu", delta = 1e-04)
   dldm <- as.vector(attr(nd, "gradient"))
   d2ldm2 <- -dldm * dldm
   d2ldm2
   }, 
            
d2ldd2 = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "sigma", delta = 1e-04)
   dldd <- as.vector(attr(nd, "gradient"))
   d2ldd2 <- -dldd*dldd
   d2ldd2
   }, 
            
d2ldv2 = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "nu", delta = 1e-04)
   dldv <- as.vector(attr(nd, "gradient"))
   d2ldv2 <- -dldv*dldv
   d2ldv2
   }, 
            
d2ldt2 = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "tau", delta = 1e-04)
   dldt <- as.vector(attr(nd, "gradient"))
   d2ldt2 <- -dldt*dldt
   },
            
d2ldmdd = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "mu", delta = 1e-04)
   dldm <- as.vector(attr(nd, "gradient"))
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "sigma", delta = 1e-04)
   dldd <- as.vector(attr(nd, "gradient"))
   d2ldmdd <- -dldm*dldd
   d2ldmdd     
   }, 
            
d2ldmdv = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "mu", delta = 1e-04)
   dldd <- as.vector(attr(nd, "gradient"))
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "nu", delta = 1e-04)
   dldv <- as.vector(attr(nd, "gradient"))
   d2ldmdv <- -dldd*dldv
   d2ldmdv  
   }, 
            
d2ldmdt = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "mu", delta = 1e-04)
   dldm <- as.vector(attr(nd, "gradient"))
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "tau", delta = 1e-04)
   dldt <- as.vector(attr(nd, "gradient"))
   d2ldmdt<- -dldm*dldt
   d2ldmdt
   }, 
            
d2ldddv = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "sigma", delta = 1e-04)
   dldd <- as.vector(attr(nd, "gradient")) 
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "nu", delta = 1e-04)
   dldv <- as.vector(attr(nd, "gradient")) 
   d2ldddv <- -dldd*dldv
   d2ldddv
   }, 
            
d2ldddt = function(y, mu, sigma, nu, tau) {
   nd = gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "sigma", delta = 1e-04)
   dldd <- as.vector(attr(nd, "gradient")) 
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "tau", delta = 1e-04)
   dldt <- as.vector(attr(nd, "gradient")) 
   d2ldddt <- -dldd*dldt
   d2ldddt
            }, 
            
d2ldvdt = function(y, mu, sigma, nu, tau) {
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "nu", delta = 1e-04)
   dldv <- as.vector(attr(nd, "gradient"))
   nd <- gamlss::numeric.deriv(dGMW(y, mu, sigma, nu, tau,
                                               log = TRUE), "tau", delta = 1e-04)
   dldt <- as.vector(attr(nd, "gradient")) 
   d2ldvdt <- -dldv*dldt
   d2ldvdt
   }, 
            
# Fin de las derivadas -------------------

   G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dGMW(y, mu, sigma, nu, tau, log = TRUE), 
        rqres = expression(rqres(pfun = "pGMW", type = "Continuous",  y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)), 
                 
   mu.initial = expression( mu <-  rep(1, length(y)) ), 
sigma.initial = expression( sigma <- rep(1, length(y)) ), 
   nu.initial = expression( nu <- rep(1, length(y)) ), 
  tau.initial = expression( tau <- rep(1, length(y)) ),
                 
     mu.valid = function(mu) all(mu > 0), 
  sigma.valid = function(sigma) all(sigma > 0), 
     nu.valid = function(nu) all(nu >= 0), 
    tau.valid = function(tau) all(tau >= 0), 
                 
      y.valid = function(y) all(y > 0)
  ), 
  class = c("gamlss.family", "family"))
}



