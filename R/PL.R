#' The Power Lindley family
#' 
#' @author Amylkar Urrea Montoya, \email{amylkar.urrea@@udea.edu.co}
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
#' @returns Returns a gamlss.family object which can be used to fit a PL distribution in the \code{gamlss()} function.
#' 
#' @example examples/examples_PL.R
#' 
#' @references
#' Almalki, S. J., & Nadarajah, S. (2014). Modifications of the 
#' Weibull distribution: A review. Reliability Engineering & 
#' System Safety, 124, 32-55.
#' 
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N., & 
#' Al-Enezi, L. J. (2013). Power Lindley distribution 
#' and associated inference. Computational Statistics & Data 
#' Analysis, 64, 20-33.
#' 
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
