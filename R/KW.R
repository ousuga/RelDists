#' The Kumaraswamy Weibull family
#' 
#' @description 
#' The Kumaraswamy Weibull distribution
#' 
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau parameter. 
#' 
#' @seealso \link{dKW}
#' 
#' @details 
#' The Kumaraswamy Weibull distribution with parameters \code{mu}, 
#' \code{sigma}, \code{nu} and \code{tau} has density given by
#' 
#' \eqn{f(x)=\nu \tau \mu \sigma x^{\sigma - 1} \exp^{-\mu x ^\sigma} [1 - \exp^{-\mu x ^ \sigma}]^{\nu-1}[1-(1-\exp^{-\mu x ^ \sigma})^\nu]^{\tau-1},}
#' 
#' for x > 0. 
#' 
#' @examples 
#' # Example 1
#' # Generating some random values with
#' # known mu, sigma, nu and tau
#' y <- rKW(n=100, mu=3, sigma=0.8, nu=2.0, tau=1.5)
#' 
#' # Fitting the model
#' require(gamlss)
#' 
#' mod <- gamlss(y~1, sigma.fo=~1, nu.fo=~1, tau.fo=~1, family='KW',
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
#' mu <- exp(2.6 + -3 * x1)
#' sigma <- exp(0.78 - 2 * x2)
#' nu <- 2
#' tau <- 1.5
#' x <- rKW(n=n, mu, sigma, nu, tau)
#' 
#' mod <- gamlss(x~x1, sigma.fo=~x2, nu.fo=~1, tau.fo=~1, family=KW,
#'               control=gamlss.control(n.cyc=5000, trace=FALSE))
#' 
#' coef(mod, what="mu")
#' coef(mod, what="sigma")
#' exp(coef(mod, what="nu"))
#' exp(coef(mod, what="tau"))
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
KW <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log") {
  mstats <- checklink("mu.link", "Kumaraswamy Weibull", 
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Kumaraswamy Weibull",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Kumaraswamy Weibull", 
                      substitute(nu.link), c("log", "own"))
  tstats <- checklink("tau.link", "Kumaraswamy Weibull", 
                      substitute(tau.link), c("log", "own"))
  
  structure(list(family = c("KW", "Kumaraswamy Weibull"), 
                 parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                 nopar = 4, 
                 type = "Continuous", 
           
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
                          
                 dldm = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/mu - y^sigma + (nu-1) * y^sigma * A / (1-A)
                   part2 <- - (tau-1) * nu * y^sigma * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldm  <- part1 + part2
                   dldm
                   },
                               
                 dldd = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/sigma + log(y) - mu*y^sigma*log(y) + (nu-1)*A*mu*y^sigma*log(y)/(1-A)
                   part2 <- -(tau-1) * mu * nu * y^sigma * log(y) * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldd  <- part1 + part2
                   dldd
                   },
                               
                 dldv = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/nu + log(1-A)
                   part2 <- -(tau-1) * (1-A)^nu * log(1-A) / (1-(1-A)^nu)
                   dldv  <- part1 + part2
                   dldv
                   },
                               
                 dldt = function(y, mu, sigma, nu, tau) {
                   A    <- exp(-mu*y^sigma)
                   dldt <- 1/tau + log(1-(1-A)^nu)
                   dldt
                   },
                               
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/mu - y^sigma + (nu-1) * y^sigma * A / (1-A)
                   part2 <- - (tau-1) * nu * y^sigma * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldm  <- part1 + part2
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                   },
                               
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/mu - y^sigma + (nu-1) * y^sigma * A / (1-A)
                   part2 <- - (tau-1) * nu * y^sigma * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldm  <- part1 + part2
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/sigma + log(y) - mu*y^sigma*log(y) + (nu-1)*A*mu*y^sigma*log(y)/(1-A)
                   part2 <- -(tau-1) * mu * nu * y^sigma * log(y) * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldd  <- part1 + part2
                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                   },
                               
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/mu - y^sigma + (nu-1) * y^sigma * A / (1-A)
                   part2 <- - (tau-1) * nu * y^sigma * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldm  <- part1 + part2
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/nu + log(1-A)
                   part2 <- -(tau-1) * (1-A)^nu * log(1-A) / (1-(1-A)^nu)
                   dldv  <- part1 + part2
                   d2ldmdv <- -dldm * dldv
                   d2ldmdv
                   },
                               
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/mu - y^sigma + (nu-1) * y^sigma * A / (1-A)
                   part2 <- - (tau-1) * nu * y^sigma * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldm  <- part1 + part2
                   A    <- exp(-mu*y^sigma)
                   dldt <- 1/tau + log(1-(1-A)^nu)
                   d2ldmdt <- -dldm * dldt
                   d2ldmdt
                   },
                               
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/sigma + log(y) - mu*y^sigma*log(y) + (nu-1)*A*mu*y^sigma*log(y)/(1-A)
                   part2 <- -(tau-1) * mu * nu * y^sigma * log(y) * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldd  <- part1 + part2
                   d2ldd2 <- -dldd * dldd
                   d2ldd2
                   },
                               
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/sigma + log(y) - mu*y^sigma*log(y) + (nu-1)*A*mu*y^sigma*log(y)/(1-A)
                   part2 <- -(tau-1) * mu * nu * y^sigma * log(y) * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldd  <- part1 + part2
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/nu + log(1-A)
                   part2 <- -(tau-1) * (1-A)^nu * log(1-A) / (1-(1-A)^nu)
                   dldv  <- part1 + part2
                   d2ldddv <- -dldd * dldv
                   d2ldddv
                   },
                               
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/sigma + log(y) - mu*y^sigma*log(y) + (nu-1)*A*mu*y^sigma*log(y)/(1-A)
                   part2 <- -(tau-1) * mu * nu * y^sigma * log(y) * A * (1-A)^(nu-1) / (1-(1-A)^nu)
                   dldd  <- part1 + part2
                   A    <- exp(-mu*y^sigma)
                   dldt <- 1/tau + log(1-(1-A)^nu)
                   d2ldddt <- -dldd * dldt
                   d2ldddt
                   },
                               
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/nu + log(1-A)
                   part2 <- -(tau-1) * (1-A)^nu * log(1-A) / (1-(1-A)^nu)
                   dldv  <- part1 + part2
                   d2ldv2 <- -dldv * dldv
                   d2ldv2
                   },
                               
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   A     <- exp(-mu*y^sigma)
                   part1 <- 1/nu + log(1-A)
                   part2 <- -(tau-1) * (1-A)^nu * log(1-A) / (1-(1-A)^nu)
                   dldv  <- part1 + part2
                   A    <- exp(-mu*y^sigma)
                   dldt <- 1/tau + log(1-(1-A)^nu)
                   d2ldvdt <- -dldv * dldt
                   d2ldvdt
                   },
                               
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   A    <- exp(-mu*y^sigma)
                   dldt <- 1/tau + log(1-(1-A)^nu)
                   d2ldt2 <- -dldt * dldt
                   d2ldt2
                   },
                             
                    G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2*dKW(y, mu, sigma, nu, tau, log=TRUE), 
                         rqres = expression(rqres(pfun="pKW", type="Continuous", 
                                               y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
                    mu.initial = expression(mu    <- rep(1, length(y))), 
                 sigma.initial = expression(sigma <- rep(1, length(y))), 
                    nu.initial = expression(nu    <- rep(1, length(y))),
                   tau.initial = expression(tau   <- rep(1, length(y))), 
                             
                      mu.valid = function(mu) all(mu > 0), 
                   sigma.valid = function(sigma) all(sigma > 0), 
                      nu.valid = function(nu) all(nu > 0), 
                     tau.valid = function(tau) all(tau > 0), 
                             
                       y.valid = function(y) all(y > 0)
            ), 
  class=c("gamlss.family", "family"))
}
