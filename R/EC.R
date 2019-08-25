#' The Exponentiated Chen family
#' 
#' @description 
#' The Exponentiated Chen family
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' 
#' @seealso \link{dEC}
#' 
#' @details 
#' The Exponentiated Chen distribution with parameters \code{mu}, 
#' \code{sigma} and \code{nu} has density given by
#' 
#' \eqn{f(x)=\mu \sigma \nu x^{\sigma - 1} \exp({x ^\sigma}) \exp\{{\nu (1 -\exp({x^\sigma})}\} 
#' \{1 - \exp({\nu (1 -\exp({x^\sigma}))})\}^{\mu-1},}
#' 
#' for \eqn{x > 0}, \eqn{\mu> 0}, \eqn{\sigma> 0}, \eqn{\nu> 0}.
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma and nu
#' y <- rEC(n=500, mu=35, sigma=0.25, nu=2)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, family='EC',
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
#' n <- 2000
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(0.044 + 3 * x1)
#' sigma <- exp(0.397 - 2 * x2)
#' nu <- 2
#' x <- rEC(n=n, mu, sigma, nu)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, family=EC,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' 
#' 
#' @references
#' \insertRef{sanku2017}{RelDists}
#'
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
EC <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Exponentiated Chen", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Chen",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Exponentiated Chen", 
                      substitute(nu.link), c("log", "own"))
  
  structure(list(family=c("EC", "Exponentiated Chen"), 
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
                   dldm <- 1 / mu + log(1 - exp(nu * (1 - exp(y^sigma))))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   b <- y^sigma
                   A <- 1 / sigma + log(y) + b * log(y) - nu * exp(b) * b * log(y)
                   B <- (mu - 1) * exp(nu * (1 - exp(b))) * nu * exp(b) * b * log(y)
                   dldd <- A + B / (1 - exp(nu * (1 - exp(b))))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   b <- y^sigma
                   C <- (1 - exp(nu * (1 - exp(b))))
                   D <- 1 / nu + (1 - exp(b)) - (mu - 1) * exp(nu * (1 - exp(b))) * (1 - exp(b)) / C
                   dldv <- D
                   dldv
                 },
                 
                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   dldm <- 1 / mu + log(1 - exp(nu * (1 - exp(y^sigma))))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   dldm <- 1 / mu + log(1 - exp(nu * (1 - exp(y^sigma))))
                   b <- y^sigma
                   A <- 1 / sigma + log(y) + b * log(y) - nu * exp(b) * b * log(y)
                   B <- (mu - 1) * exp(nu * (1 - exp(b))) * nu * exp(b) * b * log(y)
                   dldd <- A + B / (1 - exp(nu * (1 - exp(b))))
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   dldm <- 1 / mu + log(1 - exp(nu * (1 - exp(y^sigma))))
                   b <- y^sigma
                   C <- (1 - exp(nu * (1 - exp(b))))
                   D <- 1 / nu + (1 - exp(b)) - (mu - 1) * exp(nu * (1 - exp(b))) * (1 - exp(b)) / C
                   dldv <- D
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldd2 = function(y, mu, sigma, nu) {
                   b <- y^sigma
                   A <- 1 / sigma + log(y) + b * log(y) - nu * exp(b) * b * log(y)
                   B <- (mu - 1) * exp(nu * (1 - exp(b))) * nu * exp(b) * b * log(y)
                   dldd <- A + B / (1 - exp(nu * (1 - exp(b))))
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   b <- y^sigma
                   A <- 1 / sigma + log(y) + b * log(y) - nu * exp(b) * b * log(y)
                   B <- (mu - 1) * exp(nu * (1 - exp(b))) * nu * exp(b) * b * log(y)
                   dldd <- A + B / (1 - exp(nu * (1 - exp(b))))
                   C <- (1 - exp(nu * (1 - exp(b))))
                   D <- 1 / nu + (1 - exp(b)) - (mu - 1) * exp(nu * (1 - exp(b))) * (1 - exp(b)) / C
                   dldv <- D
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   b <- y^sigma
                   C <- (1 - exp(nu * (1 - exp(b))))
                   D <- 1 / nu + (1 - exp(b)) - (mu - 1) * exp(nu * (1 - exp(b))) * (1 - exp(b)) / C
                   dldv <- D
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                 },
                 
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dEC(y, mu, sigma, nu, log=TRUE), 
                 rqres = expression(rqres(pfun="pEC", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)), 
                 
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
