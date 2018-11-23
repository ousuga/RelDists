#' @name PL
#' 
#' @title 
#' The Power Lindley Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the power Lindley distribution with
#' parameters \code{mu} and \code{sigma}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu shape parameter.
#' @param sigma scale parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The power Lindley distribution with parameters \code{mu} and
#' \code{sigma} has density given by
#' 
#' f(x) = ((mu*sigma^2)/(sigma+1))*(1+x^(mu))*x^(mu-1)*exp(-sigma*(x^mu))
#' 
#' for x > 0.
#' 
#' @return 
#' \code{dPL} gives the density, \code{pPL} gives the distribution 
#' function, \code{qPL} gives the quantile function, \code{rPL}
#' generates random deviates and \code{hPL} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dPL(x, mu = 1, sigma = 0.5), from = 0, to = 15, ylim = c(0, 0.25), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pPL(x, mu = 1, sigma = 0.5), from = 0, to = 15,  ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pPL(x, mu = 1, sigma = 0.5, lower.tail = FALSE), from = 0, to = 15,  ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x=qPL(p=p, mu = 1, sigma = 0.5), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pPL(x, mu = 1, sigma = 0.5), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rPL(n = 1000, mu = 1, sigma = 0.5), freq = FALSE, , ylim = c(0, 0.25), xlab = "x", las = 1, main = "")
#' curve(dPL(x, mu = 1, sigma = 0.5),  from = 0, add = T, col = "red", ylim = c(0, 0.25))
#' 
#' ## The Hazard function
#' curve(hPL(x, mu = 1, sigma = 0.5), from = 0, to = 10, ylim = c(0, 0.5), col = "red", las = 1, ylab = "The Hazard function")
PL <- function (mu.link="log", sigma.link="log") 
{
  mstats <- checklink("mu.link", "Power Lindley", substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Power Lindley", substitute(sigma.link), c("log", "own"))
  
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
                 
                 dldm = function(y,mu,sigma) ((1/mu)+(((y^mu)*log(y))/(1+(y^mu)))+log(y)-(sigma*(y^mu)*log(y))),
                 d2ldm2 = function(y,mu,sigma) {
                   dldm = function(y,mu,sigma) (-(1/(mu^2))+ ((y^mu)*(log(y))^2)/(1+(y^mu)^2) -(sigma*(y^mu)*(log(y))^2) ) 
                   ans <- dldm(y,mu,sigma)
                   ans <- -ans^2
                 },
                 dldd = function(y,mu,sigma) ((2/sigma)-(1/(sigma+1))-(y^mu)),
                 d2ldd2 = function(y,mu,sigma) {
                   dldd = function(y,mu,sigma) ((-2/(sigma^2))-((1/(sigma+1))^2))
                   ans <- dldd(y,mu,sigma)
                   ans <- -ans^2
                 },
                 
                 d2ldmdd = function(y,mu,sigma) -(-(y^mu)*log(y))^2,
                 
                 G.dev.incr  = function(y,mu,sigma,...) -2*dPL(y, mu, sigma, log=TRUE), 
                 rqres = expression(rqres(pfun="pPL", type="Continuous", y=y, mu=mu, sigma=sigma)),
                 
                 mu.initial = expression( mu <-  rep(1, length(y)) ),     
                 sigma.initial = expression( sigma <- rep(1, length(y)) ), 
                 
                 mu.valid = function(mu) all(mu > 0) , 
                 sigma.valid = function(sigma)  all(sigma > 0), 
                 
                 y.valid = function(y)  all(y > 0)
  ),
  class = c("gamlss.family","family"))
}
#' @export
#' @rdname PL
dPL<-function(x,mu,sigma, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  loglik<- log(mu) + 2*log(sigma) - log(sigma+1) +
    log(1+(x^mu)) + (mu-1)*log(x) - sigma*(x^mu)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname PL
pPL <- function(q,mu,sigma, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  cdf <- 1-(1+((sigma/(sigma+1))*q^mu))*exp(-sigma*(q^mu))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname PL
qPL <- function(p,mu,sigma, lower.tail = TRUE, log.p = FALSE){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  fda <- function(x,mu, sigma){
    1-(1+((sigma/(sigma+1))*x^mu))*exp(-sigma*(x^mu))
  }
  
  fda1 <- function(x, mu, sigma, p) {fda(x, mu, sigma) - p}
  
  r_de_la_funcion <- function(mu, sigma, p) {
    uniroot(fda1, interval=c(0,1e+06), mu, sigma, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, p)
  q
}


#' @export
#' @rdname PL
rPL <- function(n,mu,sigma){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qPL(p, mu,sigma)
  r
}

#' @export
#' @rdname PL
hPL<-function(x,mu,sigma){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma<=0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  
  h <- dPL(x,mu,sigma, log = FALSE)/pPL(q=x,mu,sigma, lower.tail=FALSE, log.p = FALSE)
  h  
}
