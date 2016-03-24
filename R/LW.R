#' @name LW
#' 
#' @title 
#' The Log-Weibull Distribution
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the log-weibull distribution with
#' parameters \code{a} and \code{b}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param a parameter one.
#' @param b parameter two.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#'  
#' @details 
#' The log-weibull distribution with parameters \code{a} and 
#' \code{b} has density given by
#' 
#' f(x)= (1/b)*exp((x-a)/b)*exp(-exp((x-a)/b))
#' 
#' for -Inf < x < Inf.
#' 
#' @return 
#' \code{dLW} gives the density, \code{pLW} gives the distribution 
#' function, \code{qLW} gives the quantile function, \code{rLW}
#' generates random deviates and \code{hLW} gives the hazard function.
#' 
#' @export
#' @examples  
#' ## The probability density function 
#' curve(dLW(x, a = 0, b = 1), from = -20, to = 10, ylim = c(0, 0.4), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pLW(x, a = 0, b = 1), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pLW(x, a = 0, b = 1, lower.tail = FALSE), from = 0, to = 10, ylim = c(0, 1), col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function
#' p <- seq(from = 0, to = 0.998, length.out = 100)
#' plot(x=qLW(p, a = 0, b = 1), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pLW(x, a = 0, b = 1), from = -6, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rLW(10000, a = 0, b = 1), freq = FALSE,  ylim = c(0, 0.4), xlab = "x", las = 1, main = "")
#' curve(dLW(x, a = 0, b = 1),  from = -20, to = 10,  ylim = c(0, 0.4), add = TRUE, col = "red") 
#' 
#' ## The Hazard function
#' curve(hLW(x, a = 0, b = 1), from = -20, to = 0, ylim = c(0, 0.3), col = "red", ylab = "The hazard function", las = 1)
#' 
dLW<-function(x,a,b, log = FALSE){
  if (any(b<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  loglik<- -log(b) + (x-a)/b - exp((x-a)/b)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname LW
pLW <- function(q,a,b, lower.tail=TRUE, log.p = FALSE){
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  cdf <- 1-exp(-exp((q-a)/b))
  
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname LW
qLW <- function(p,a,b, lower.tail = TRUE, log.p = FALSE){
  
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- b*(log(-log(1-p)))+a  
  q
}
#' @export
#' @rdname LW
rLW <- function(n,a,b){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qLW(p,a,b)
  r
}
#' @export
#' @rdname LW
hLW<-function(x,a,b){
  if (any(b<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  
  h <- dLW(x,a,b, log = FALSE)/pLW(q=x,a,b, lower.tail=FALSE, log.p = FALSE)
  h
}
