#' The Muth family
#' 
#' @description 
#' The Muth family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' 
#' @seealso \link{dMuth}
#' 
#' @details 
#' The Muth distribution with parameter \code{mu} has density given by
#' 
#' \eqn{f(x)= ({exp(\mu x) - \mu})exp[{\mu x - \frac{1}{\mu} ({exp(\mu x) - 1})}]}
#' 
#' for \eqn{x > 0}, \eqn{0< \mu <=1}.
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu
#' y <- rMuth(n=200, mu=0.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, family='Muth',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu
#' # using the inverse link function
#' exp(coef(mod, what='mu'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 2000
#' x1 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(-2.19 + 3 * x1)
#' x <- rMuth(n=n, mu)
#' 
#' mod <- gamlss(x~x1, family=Muth,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' 
#' @references
#' \insertRef{abdullah2018}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' 
#' @export
Muth <- function (mu.link="log") {
  mstats <- checklink("mu.link", "Muth", 
                      substitute(mu.link), c("logit", "own"))
  
  structure(list(family=c("Muth", "Muth"), 
                 parameters=list(mu=TRUE), 
                 nopar=1, 
                 type="Continuous", 
                 
                 mu.link = as.character(substitute(mu.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 
                 mu.dr = mstats$mu.eta, 
                 
                 # Primeras derivadas ---------------------------------
                 dldm = function(y, mu) {
                   t <- exp(mu * y)
                   A <- (t * y - 1) / (t - mu)
                   dldm <- A + y - (mu * t * y - t + 1) / mu^2
                   dldm
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu) {
                   t <- exp(mu * y)
                   A <- (t * y - 1) / (t - mu)
                   dldm <- A + y - (mu * t * y - t + 1) / mu^2
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 
                 G.dev.incr = function(y, mu, ...) -2*dMuth(y, mu, log=TRUE), 
                 rqres = expression(rqres(pfun="pMuth", type="Continuous", y=y, mu=mu)), 
                 
                 mu.initial = expression(mu       <- rep(1, length(y))), 
                 
                 mu.valid = function(mu)       all(mu > 0 & mu <= 1),
                 
                 y.valid = function(y) all(y > 0)
  ), 
  class=c("gamlss.family", "family"))
}
