#' The Extended Weibull family
#' 
#' @description 
#' The Extended Weibull family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dExW}
#' 
#' @details 
#' The Extended Weibull distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma \nu x^{\sigma -1} exp({-\mu x^{\sigma}})} {[1 -(1-\nu) exp({-\mu x^{\sigma}})]^2},}
#' 
#' for x > 0. 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rExW(n=200, mu=0.3, sigma=2, nu=0.05)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='ExW',
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
#' n <- 500
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(-2 + 3 * x1)
#' sigma <- exp(1.3 - 2 * x2)
#' nu <- 0.05
#' x <- rExW(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=ExW,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Zhang2007}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
ExW <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Extended Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Extended Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Extended Weibull", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("ExW", "Extended Weibull"), 
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
                    A    <- y^(sigma) * exp(- mu * y^sigma)
                    B    <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldm <- 1 / mu - y^sigma + (2 * A) / B 
                    dldm
                  },
                  
                  dldd = function(y, mu, sigma, nu) {
                    A    <- 1 / sigma + log(y) - mu * y^(sigma) * log(y)
                    B    <- exp(- mu * y^sigma) * mu * y^(sigma) * log(y)
                    C    <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldd <- A + (2 * B) / C
                    dldd
                  },
                  
                  dldv = function(y, mu, sigma, nu) {
                    A    <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldv <- 1 / nu - (2 * exp(- mu * y^sigma))/ A 
                    dldv
                  },
                  
                     
                 # Segundas derivadas ---------------------------------
                  d2ldm2 = function(y, mu, sigma, nu) {
                    A      <- y^(sigma) * exp(- mu * y^sigma)
                    B      <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldm   <- 1 / mu - y^sigma + (2 * A) / B 
                    d2ldm2 <- -dldm * dldm
                    d2ldm2
                  },
                  
                  d2ldmdd = function(y, mu, sigma, nu) {
                    A       <- y^(sigma) * exp(- mu * y^sigma)
                    B       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldm    <- 1 / mu - y^sigma + (2 * A) / B 
                    C       <- 1 / sigma + log(y) - mu * y^(sigma) * log(y)
                    D       <- exp(- mu * y^sigma) * mu * y^(sigma) * log(y)
                    E       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldd    <- C + (2 * D) / E
                    d2ldmdd <- -dldm * dldd
                    d2ldmdd
                  },
                  
                  d2ldmdv = function(y, mu, sigma, nu) {
                    A       <- y^(sigma) * exp(- mu * y^sigma)
                    B       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldm    <- 1 / mu - y^sigma + (2 * A) / B 
                    C       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldv    <- 1 / nu - (2 * exp(- mu * y^sigma))/ C
                    d2ldmdv <- -dldm * dldv
                    d2ldmdv
                  },
                  
                  d2ldd2 = function(y, mu, sigma, nu) {
                    A      <- 1 / sigma + log(y) - mu * y^(sigma) * log(y)
                    B      <- exp(- mu * y^sigma) * mu * y^(sigma) * log(y)
                    C      <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldd   <- A + (2 * B) / C
                    d2ldd2 <- -dldd * dldd
                    d2ldd2
                  },
                  
                  d2ldddv = function(y, mu, sigma, nu) {
                    A       <- 1 / sigma + log(y) - mu * y^(sigma) * log(y)
                    B       <- exp(- mu * y^sigma) * mu * y^(sigma) * log(y)
                    C       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldd    <- A + (2 * B) / C
                    D       <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldv    <- 1 / nu - (2 * exp(- mu * y^sigma))/ D 
                    d2ldddv <- -dldd * dldv
                    d2ldddv
                  },
                  
                  d2ldv2 = function(y, mu, sigma, nu) {
                    A      <- (1 - (1 - nu) * exp(- mu * y^sigma))
                    dldv   <- 1 / nu - (2 * exp(- mu * y^sigma))/ A 
                    d2ldv2 <- -dldv * dldv
                    d2ldv2
                  },
              
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dExW(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pExW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
                 mu.initial = expression(mu       <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 nu.initial = expression(nu       <- rep(1, length(y))),
                 
                 mu.valid = function(mu)       all(mu > 0), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu)       all(nu > 0),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
