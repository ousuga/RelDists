#' The Ex-Wald family
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' The function \code{ExWALD()} defines the Ex-wALD distribution, three-parameter
#' continuous distribution for a \code{gamlss.family} object to be used in 
#' GAMLSS fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma parameter.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu parameter.
#'
#' @references
#' Schwarz, W. (2001). The ex-Wald distribution as a descriptive model 
#' of response times. Behavior Research Methods, 
#' Instruments, & Computers, 33, 457-469.
#' 
#' Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to 
#' response time data: An example using functions for the S-PLUS package. 
#' Behavior Research Methods, Instruments, & Computers, 36, 678-694.
#'
#' @seealso \link{dExWALD}.
#'
#' @details
#' The Ex-Wald distribution with parameters \eqn{\mu}, \eqn{\sigma} 
#' and \eqn{\nu} has density given by
#'
#' \eqn{f(x |\mu, \sigma, \nu) = \frac{1}{\nu} \exp(\frac{-x}{\nu} + \sigma(\mu-k)) F_W(x|k, \sigma) \, \text{for} \, k \geq 0}
#'
#' \eqn{f(x |\mu, \sigma, \nu) = \frac{1}{\nu} \exp\left( \frac{-(\sigma-\mu)^2}{2x} \right) Re \left( w(k^\prime \sqrt{x/2} + \frac{\sigma i}{\sqrt{2x}}) \right)  \, \text{for} \, k < 0}
#'
#' where \eqn{k=\sqrt{\mu^2-\frac{2}{\nu}}}, 
#' \eqn{k^\prime=\sqrt{\frac{2}{\nu}-\mu^2}} and
#' \eqn{F_W} corresponds to the cumulative function of 
#' the Wald distribution. 
#' 
#' More details about those expressions
#' can be found on page 680 from Heathcote (2004).
#'
#' @returns Returns a gamlss.family object which can be used to fit a 
#' Ex-WALD distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_ExWALD.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
ExWALD <- function(mu.link="log", 
                   sigma.link="log", 
                   nu.link="log") {
  
  mstats <- checklink("mu.link", "Ex-Wald",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "Ex-Wald",
                      substitute(sigma.link), c("log"))
  vstats <- checklink("nu.link", "Ex-Wald",
                      substitute(nu.link), c("log"))
  
  structure(list(family     = c("ExWALD", "Ex-Wald"),
                 parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
                 nopar      = 3,
                 type       = "Continuous",
                 
                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link    = as.character(substitute(nu.link)),
                 
                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 nu.linkfun    = vstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 nu.linkinv    = vstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 nu.dr    = vstats$mu.eta,
                 
                 # First derivates
                 
                 dldm = function(y, mu, sigma, nu) {
                   dm <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                               theta="mu",
                                               delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma, nu) {
                   dd <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                               theta="sigma",
                                               delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },
                 
                 dldv = function(y, mu, sigma, nu) {
                   dv <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                               theta="nu",
                                               delta=0.001)
                   dldv <- as.vector(attr(dv, "gradient"))
                   dldv
                 },

                 # Second derivates
                 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   dm <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm   <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm*dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma, nu) {
                   dm <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd    <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm*dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   dm <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dv   <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.001)
                   dldv    <- as.vector(attr(dv, "gradient"))
                   d2ldmdv <- - dldm*dldv
                   d2ldmdv <- ifelse(d2ldmdv < -1e-15, d2ldmdv, -1e-15)
                   d2ldmdv
                 },
                 
                 d2ldd2  = function(y, mu, sigma, nu) {
                   dd  <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd   <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd*dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, sigma, nu) {
                   dd <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dv   <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.001)
                   dldv    <- as.vector(attr(dv, "gradient"))
                   d2ldddv <- - dldd*dldv
                   d2ldddv <- ifelse(d2ldddv < -1e-15, d2ldddv, -1e-15)
                   d2ldddv
                 },
                 
                 d2ldv2  = function(y, mu, sigma, nu) {
                   dv   <- gamlss::numeric.deriv(dExWALD(y, mu, sigma, nu, log=TRUE),
                                                 theta="nu",
                                                 delta=0.001)
                   dldv   <- as.vector(attr(dv, "gradient"))
                   d2ldv2 <- - dldv*dldv
                   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
                   d2ldv2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dExWALD(y, mu, sigma, nu, log=TRUE),
                 rqres      = expression(rqres(pfun="pExWALD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
                 
                 mu.initial    = expression(mu    <- rep(fitexw(y)[1], length(y))),
                 sigma.initial = expression(sigma <- rep(fitexw(y)[2], length(y))),
                 nu.initial    = expression(nu    <- rep(fitexw(y)[3], length(y))),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 nu.valid    = function(nu)    all(nu > 0),
                 
                 y.valid = function(y) all(y > 0),
                 
                 mean     = function(mu, sigma, nu) {nu + sigma/mu}, 
                 variance = function(mu, sigma, nu) {nu^2 + sigma/mu^3}
  ),
  class = c("gamlss.family", "family"))
}
#' Auxiliar function for the Ex-Wald distribution
#' @description This function generates start values.
#' @param x a value for x.
#' @param p a value for p.
#' @return returns a vector with starting values.
#' @keywords internal
#' @importFrom stats sd var
#' @export
exwstpt <- function(x, p=0.5) {
  t <- p*sd(x)
  m <- sqrt((mean(x)-t)/(var(x)-t^2))
  a <- m*(mean(x)-t)
  res <- c(m, a, t)
  return(res)
}
#' Auxiliar pwald function for the Ex-Wald distribution
#' @description This function emulates the pWALD function.
#' @param w a value for w.
#' @param m a value for m.
#' @param a a value for a.
#' @param s a value for s.
#' @return returns a vector with cumulative probabilities.
#' @keywords internal
#' @export
pwald <- function(w, m, a, s=0) {
  w <- w - s; 
  sqrtw <- sqrt(w); 
  k1 <- (m*w-a)/sqrtw; 
  k2 <- (m*w+a)/sqrtw
  p1 <- exp(2*a*m); 
  p2 <- pnorm(-k2); 
  bad <- (p1==Inf) | (p2==0); 
  p <- p1*p2
  p[bad] <- (exp(-(k1[bad]^2)/2 - 0.94/(k2[bad]^2))/(k2[bad]*((2*pi)^.5)))
  p + pnorm(k1)
}
#' Auxiliar function to generate random values for Ex-Wald distribution
#' @description This function generates another form to generate values.
#' @param x a value for x.
#' @param y a value for y.
#' @return returns a vector with values.
#' @keywords internal
#' @export
rew <- function(x, y) {
  uv <- uandv(y, x)
  exp(y^2-x^2)*(cos(2*x*y)*(1-uv[,1])+sin(2*x*y)*uv[,2])
}
#' Auxiliar den_exw function for the Ex-Wald distribution
#' @description This function emulates the dExWALD function.
#' @param r value for r.
#' @param m value for m.
#' @param a value for a.
#' @param t value for t.
#' @return returns a vector with probabilities.
#' @keywords internal
#' @export
den_exw <- function(r, m, a, t) {
  k <- (m^2 - (2/t))
  if (k < 0) {
    res <- exp(m*a - (a^2)/(2*r) - r*(m^2)/2)*
      rew(sqrt(-r*k/2), a/sqrt(2*r))/t
  } else {
    k <- sqrt(k)
    res <- pwald(r,k,a)*exp(a*(m-k) - (r/t))/t
  }
  return(res)
}
#' Auxiliar function to obtain minus loglik for Ex-Wald 
#' @description This function calculates minus loglik.
#' @param p a vector with the parameters.
#' @param x a vector with the random sample.
#' @return returns the minus loglik.
#' @keywords internal
#' @export
negllexw <- function(p, x) {
  -sum(log(den_exw(x, p[1], p[2], p[3])))
}
#' Auxiliar function for the Ex-Wald distribution
#' @description This function generates starting values.
#' @param rt a vector with the random sample.
#' @param p a value for p.
#' @param start an optional start vector.
#' @param scaleit logical value to scale.
#' @return returns a vector with cumulative probabilities.
#' @keywords internal
#' @importFrom stats nlminb
#' @export
fitexw <- function(rt, p=0.5,
                   start=exwstpt(rt, p),
                   scaleit=TRUE){
  if (scaleit) 
    fit <- nlminb(start=start,
                  lower=c(1e-8, 1e-8, 1),
                  objective=negllexw,
                  x=rt,
                  scale=(1/start),
                  control=list(eval.max=400, iter.max=300))
  else 
    fit <- nlminb(start=start,
                  lower=c(1e-8, 1e-8, 1),
                  objective=negllexw,
                  x=rt,
                  control=list(eval.max=400, iter.max=300, 
                               scale.tol=1e-8))

  return(fit$par)
}


