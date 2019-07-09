#' The Additive Weibull family
#' 
#' @description 
#' The Additive Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dAddW}
#' 
#' @details 
#' Additive Weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x) = (\mu\nu x^{\nu - 1} + \sigma\tau x^{\tau - 1}) \exp({-\mu x^{\nu} - \sigma x^{\tau} }),}
#' 
#' for x > 0. 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rAddW(n=100, mu=1.5, sigma=0.2, nu=3, tau=0.8)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family='AddW',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu, sigma, nu and tau
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' exp(coef(mod, what='sigma'))
#' exp(coef(mod, what='nu'))
#' exp(coef(mod, what='tau'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 200
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(1.67 + -3 * x1)
#' sigma <- exp(0.69 - 2 * x2)
#' nu <- 3
#' tau <- 0.8
#' x <- rAddW(n=n, mu, sigma, nu, tau)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=AddW,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' exp(coef(mod, what="tau"))
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Xie1996}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
AddW <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log") {
  mstats <- checklink("mu.link", "Additive Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Additive Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Additive Weibull", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Additive Weibull", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family=c("AddW", "Additive Weibull"), 
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                 nopar=4, 
                 type="Continuous", 
                 
                 mu.link = as.character(substitute(mu.link)), 
              sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)),
                tau.link = as.character(substitute(tau.link)),
                 
              mu.linkfun = mstats$linkfun, 
           sigma.linkfun = dstats$linkfun, 
              nu.linkfun = vstats$linkfun,
             tau.linkfun = tstats$linkfun,
                 
              mu.linkinv = mstats$linkinv, 
           sigma.linkinv = dstats$linkinv, 
              nu.linkinv = vstats$linkinv,
             tau.linkinv = tstats$linkinv,
                 
                   mu.dr = mstats$mu.eta, 
                sigma.dr = dstats$mu.eta, 
                   nu.dr = vstats$mu.eta,
                  tau.dr = tstats$mu.eta,
                 
           # Primeras derivadas ---------------------------------
           dldm = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm <- (nu * y^(nu - 1))/A - y^nu
             dldm
           },
           
           dldd = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd <- (tau * y^(tau - 1))/A - y^tau
             dldd
           },
           
           dldv = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B    <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv <- (mu * B)/A - mu * log(y) * y^nu 
             dldv
           },
           
           dldt = function(y, mu, sigma, nu, tau) {
             A    <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B    <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt <- (sigma * B)/A - sigma * log(y) * y^tau 
             dldt
           },
           
           # Segundas derivadas ---------------------------------
           d2ldm2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm   <- (nu * y^(nu - 1))/A - y^nu
             d2ldm2 <- -dldm * dldm
             d2ldm2
           },
           
           d2ldmdd = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             d2ldmdd <- -dldm * dldd
             d2ldmdd
           },
           
           d2ldmdv = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldmdv <- -dldm * dldv
             d2ldmdv
           },
           
           d2ldmdt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldm    <- (nu * y^(nu - 1))/A - y^nu
             B       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldmdt <- -dldm * dldt
             d2ldmdt
           },
           
           d2ldd2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd   <- (tau * y^(tau - 1))/A - y^tau
             d2ldd2 <- -dldd * dldd
             d2ldd2
           },
           
           d2ldddv = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldddv <- -dldd * dldv
             d2ldddv
           },
           
           d2ldddt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             dldd    <- (tau * y^(tau - 1))/A - y^tau
             B       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldddt <- -dldd * dldt
             d2ldddt
           },
           
           d2ldv2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B      <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv   <- (mu * B)/A - mu * log(y) * y^nu 
             d2ldv2 <- -dldv * dldv
             d2ldv2
           },
           
           d2ldvdt = function(y, mu, sigma, nu, tau) {
             A       <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B       <- y^(nu - 1) + nu * log(y) * y^(nu - 1)
             dldv    <- (mu * B)/A - mu * log(y) * y^nu 
             C       <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt    <- (sigma * C)/A - sigma * log(y) * y^tau 
             d2ldvdt <- -dldv * dldt
             d2ldvdt
           },
           
           d2ldt2 = function(y, mu, sigma, nu, tau) {
             A      <- mu * nu * y^(nu - 1) + sigma * tau * y^(tau - 1)
             B      <- y^(tau - 1) + tau * log(y) * y^(tau - 1)
             dldt   <- (sigma * B)/A - sigma * log(y) * y^tau 
             d2ldt2 <- -dldt * dldt
             d2ldt2
           },
           
            G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dAddW(y, mu, sigma, nu, tau, log=TRUE), 
                 rqres = expression(rqres(pfun="pAddW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                 
            mu.initial = expression(mu    <- rep(1, length(y))), 
         sigma.initial = expression(sigma <- rep(1, length(y))), 
            nu.initial = expression(nu    <- rep(1, length(y))),
           tau.initial = expression(tau    <- rep(1, length(y))),
                 
              mu.valid = function(mu)    all(mu > 0), 
           sigma.valid = function(sigma) all(sigma > 0), 
              nu.valid = function(nu)    all(nu > 0),
             tau.valid = function(tau)    all(tau > 0),
                 
               y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}