#' The Log-Weibull family
#' 
#' @description 
#' The Log-Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @seealso \link{dLW}
#' 
#' @details 
#' The Log-Weibull Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(y)=(1/\sigma) e^{((y - \mu)/\sigma)} exp\{-e^{((y - \mu)/\sigma)}\},}
#' 
#' for - \code{infty} < y < \code{infty}.
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu and sigma 
#' y <- rLW(n=100, mu=0, sigma=1.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, family= 'LW',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu and sigma
#' # using the inverse link function
#' coef(mod, 'mu')
#' exp(coef(mod, 'sigma'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 200
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- 1.5 - 3 * x1
#' sigma <- exp(1.4 - 2 * x2)
#' x <- rLW(n=n, mu, sigma)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, family=LW,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Gumbel1958}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
LW <- function (mu.link="identity", sigma.link="log") {
  mstats <- checklink("mu.link", "Log-Weibull", 
                      substitute(mu.link), c("identity", "own"))
  dstats <- checklink("sigma.link", "Log-Weibull",
                      substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("LW", "Log-Weibull"), 
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
                 
                 dldm = function(y, mu, sigma) {
                   dldm <- (-1)/sigma + exp((y-mu)/sigma)/sigma
                   dldm
                   },
                               
                 dldd = function(y, mu, sigma) {
                   A     <- exp((y-mu)/sigma)
                   part1 <- (-1)/sigma - (y-mu)/(sigma * sigma)
                   part2 <- ((y-mu)/(sigma * sigma)) * A
                   dldd  <- part1 + part2
                   dldd
                   },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   dldm   <- (-1)/sigma + exp((y-mu)/sigma)/sigma
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                   },
                               
                 d2ldmdd = function(y, mu, sigma) {
                   dldm    <- (-1)/sigma + exp((y-mu)/sigma)/sigma
                   A       <- exp((y-mu)/sigma)
                   part1   <- (-1)/sigma - (y-mu)/(sigma * sigma)
                   part2   <- ((y-mu)/(sigma * sigma)) * A
                   dldd    <- part1 + part2
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                   },
                               
                 d2ldd2 = function(y, mu, sigma) {
                   A      <- exp((y-mu)/sigma)
                   part1  <- (-1)/sigma - (y-mu)/(sigma * sigma)
                   part2  <- ((y-mu)/(sigma * sigma)) * A
                   dldd   <- part1 + part2
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                   },
                               
                    G.dev.incr = function(y, mu, sigma, ...) -2*dLW(y, mu, sigma, log=TRUE), 
                         rqres = expression(rqres(pfun="pLW", type="Continuous", y=y, mu=mu, sigma=sigma)), 
                             
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                             
                      mu.valid = function(mu) TRUE, 
                   sigma.valid = function(sigma) all(sigma > 0), 
                             
                       y.valid = function(y) TRUE
    ), 
  class=c("gamlss.family", "family"))
}
