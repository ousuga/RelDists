#' The Power Lindley family
#' 
#' @description 
#' Power Lindley distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' 
#' @seealso \link{dPL}
#' 
#' @details 
#' The Power Lindley Distribution with parameters \code{mu} 
#' and \code{sigma} has density given by
#' 
#' \eqn{f(x) = \frac{\mu \sigma^2}{\sigma + 1} (1 + x^\mu) x ^ {\mu - 1} \exp({-\sigma x ^\mu}),}
#' 
#' for x > 0.
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu and sigma 
#' y <- rPL(n=100, mu=1.5, sigma=0.2)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, family= 'PL',
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' # Extracting the fitted values for mu and sigma
#' # using the inverse link function
#' exp(coef(mod, 'mu'))
#' exp(coef(mod, 'sigma'))
#' 
#' # Example 2
#' # Generating random values under some model
#' n <- 200
#' x1 <- runif(n, min=0.4, max=0.6)
#' x2 <- runif(n, min=0.4, max=0.6)
#' mu <- exp(1.2 - 2 * x1)
#' sigma <- exp(0.8 - 3 * x2)
#' x <- rPL(n=n, mu, sigma)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, family=PL,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' 
#' @references
#' \insertRef{almalki2014modifications}{RelDists}
#' 
#' \insertRef{Ghitanya2013}{RelDists}
#' 
#' @importFrom Rdpack reprompt
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
PL <- function (mu.link="log", sigma.link="log") {
  mstats <- checklink("mu.link", "Power Lindley", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Power Lindley",
                      substitute(sigma.link), c("log", "own"))
  
  structure(list(family = c("PL", "Power Lindley"), 
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
                   part1 <- log(y) - sigma * y^mu * log(y)
                   part2 <- 1/mu + (y^(mu) * log(y))/(1+y^mu) 
                   dldm  <- part1 + part2
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   dldd <- 2/sigma - 1/(sigma+1) - y^mu
                   dldd
                 },
                 
                 d2ldm2 = function(y, mu, sigma) {
                   part1  <- log(y) - sigma * y^mu * log(y)
                   part2  <- 1/mu + (y^(mu) * log(y))/(1+y^mu)
                   dldm   <- part1 + part2
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   part1   <- log(y) - sigma * y^mu * log(y)
                   part2   <- 1/mu + (y^(mu) * log(y))/(1+y^mu) 
                   dldm    <- part1 + part2
                   dldd    <- 2/sigma - 1/(sigma+1) - y^mu
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldd2 = function(y, mu, sigma) {
                   dldd   <- 2/sigma - 1/(sigma+1) - y^mu
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                 },
                 
                    G.dev.incr = function(y, mu, sigma, ...) -2*dPL(y, mu, sigma, log=TRUE), 
                         rqres = expression(rqres(pfun="pPL", type="Continuous", y=y, mu=mu, sigma=sigma)), 
                 
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                 
                      mu.valid = function(mu) all(mu > 0), 
                   sigma.valid = function(sigma) all(sigma > 0), 
                 
                       y.valid = function(y) all(y > 0)
    ), 
    class=c("gamlss.family", "family"))
}
