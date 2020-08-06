#' The Generalized Modified Weibull family 
#' 
#' @description 
#' The Generalized modified Weibull distribution
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
#' @example examples/examples_GMW.R
#' 
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GMW <- function (mu.link = "log", sigma.link = "log", nu.link = "sqrt", tau.link = "sqrt") {
  mstats <- checklink("mu.link",    "Generalized Modified Weibull", substitute(mu.link),    c("log", "own"))
  dstats <- checklink("sigma.link", "Generalized Modified Weibull", substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link",    "Generalized Modified Weibull", substitute(nu.link),    c("sqrt", "own"))
  tstats <- checklink("tau.link",   "Generalized Modified Weibull", substitute(tau.link),    c("sqrt", "own"))
  
  structure(list(family = c("GMW", "Generalized Modified Weibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE), 
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
                   dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/
                                                       (1-exp(-mu*(y^nu)*exp(tau*y))))
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
                 
                 d2ldm2 = function(y, mu, sigma, nu, tau) {
                   dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/
                                                       (1-exp(-mu*(y^nu)*exp(tau*y))))
                    d2ldm2 <- -dldm * dldm
                    d2ldm2 
                    }, 
                          
                 d2ldmdd = function(y, mu, sigma, nu, tau) {
                    dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/
                                                         (1-exp(-mu*(y^nu)*exp(tau*y))))
                    dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
                    d2ldmdd <- -dldm*dldd
                    d2ldmdd     
                    }, 
                             
                 d2ldmdv = function(y, mu, sigma, nu, tau) {
                   dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
                   dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/
                                                       (1-exp(-mu*(y^nu)*exp(tau*y))))
                   part1 <- mu*exp(tau*y)*(y^nu)*log(y)
                   part2 <- mu*(y^nu)*exp(tau*y)
                   dldv  <- log(y) + 1/(nu + tau*y) - part1 + (((sigma-1)*part1*exp(-part2))/(1-exp(-part2)))
                   d2ldmdv <- -dldd*dldv
                   d2ldmdv  
                   }, 
            
                 d2ldmdt = function(y, mu, sigma, nu, tau) {
                   dldm <- (1/mu)-(y^nu)*exp(tau*y)+(((sigma-1)*(y^nu)*exp(tau*y)*exp(-mu*(y^nu)*exp(tau*y)))/
                                                       (1-exp(-mu*(y^nu)*exp(tau*y))))
                   part1 <- mu*(y^nu)*exp(tau*y)
                   dldt  <- (y/(nu+(tau*y))) + y -y*part1 +(((sigma-1)*exp(-part1)*y*part1)/(1-exp(-part1)))
                   d2ldmdt <- -dldm*dldt
                   d2ldmdt
                   }, 
                 
                 d2ldd2 = function(y, mu, sigma, nu, tau) {
                   dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
                   d2ldd2 <- -dldd*dldd
                   d2ldd2
                   }, 
                 
                 d2ldddv = function(y, mu, sigma, nu, tau) {
                   dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
                   part1 <- mu*exp(tau*y)*(y^nu)*log(y)
                   part2 <- mu*(y^nu)*exp(tau*y)
                   dldv  <- log(y) + 1/(nu + tau*y) - part1 + (((sigma-1)*part1*exp(-part2))/(1-exp(-part2)))
                   d2ldddv <- -dldd*dldv
                   d2ldddv
                   }, 
            
                 d2ldddt = function(y, mu, sigma, nu, tau) {
                   dldd <- (1/sigma)+log(1-exp(-mu*(y^nu)*exp(tau*y)))
                   part1 <- mu*(y^nu)*exp(tau*y)
                   dldt  <- (y/(nu+(tau*y))) + y -y*part1 +(((sigma-1)*exp(-part1)*y*part1)/(1-exp(-part1)))
                   d2ldddt <- -dldd*dldt
                   d2ldddt
                   }, 
                 
                 d2ldv2 = function(y, mu, sigma, nu, tau) {
                   part1 <- mu*exp(tau*y)*(y^nu)*log(y)
                   part2 <- mu*(y^nu)*exp(tau*y)
                   dldv  <- log(y) + 1/(nu + tau*y) - part1 + (((sigma-1)*part1*exp(-part2))/(1-exp(-part2)))
                   d2ldv2 <- -dldv*dldv
                   d2ldv2
                   }, 
            
                 d2ldvdt = function(y, mu, sigma, nu, tau) {
                   part1 <- mu*exp(tau*y)*(y^nu)*log(y)
                   part2 <- mu*(y^nu)*exp(tau*y)
                   dldv  <- log(y) + 1/(nu + tau*y) - part1 + (((sigma-1)*part1*exp(-part2))/(1-exp(-part2)))
                   part1 <- mu*(y^nu)*exp(tau*y)
                   dldt  <- (y/(nu+(tau*y))) + y -y*part1 +(((sigma-1)*exp(-part1)*y*part1)/(1-exp(-part1)))
                   d2ldvdt <- -dldv*dldt
                   d2ldvdt
                   }, 
                 
                 d2ldt2 = function(y, mu, sigma, nu, tau) {
                   part1 <- mu*(y^nu)*exp(tau*y)
                   dldt  <- (y/(nu+(tau*y))) + y -y*part1 +(((sigma-1)*exp(-part1)*y*part1)/(1-exp(-part1)))
                   d2ldt2 <- -dldt*dldt
                   },
                 
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




