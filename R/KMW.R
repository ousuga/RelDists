#' @name KMW
#' 
#' @title 
#' The Kumaraswamy modified Weibull distribution 
#' 
#' @description 
#' Density, distribution function, quantile function, 
#' random generation and hazard function for the kumaraswamy modified weibull distribution with
#' parameters \code{alpha}, \code{theta}, \code{lambda}, \code{a} and \code{b}.
#' 
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param alpha parameter one.    
#' @param beta parameter two.
#' @param lambda parameter three.
#' @param a parameter four.
#' @param b parameter five.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).	
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' P[X <= x], otherwise, P[X > x].
#' 
#' @details 
#' The kumaraswamy modified weibull with parameters \code{alpha}, \code{theta}, 
#' \code{lambda}, \code{a} and \code{b} has density given by
#' 
#' f(x)=a*b*alpha*x^(theta-1)*(theta+lambda*x)*exp(lambda*x-alpha*x^(theta)*exp(lambda*x))*
#' (1-exp(-alpha*x^theta*exp(lambda*x)))^(a-1)*
#' (1-(1-exp(-aplha*x^theta*exp(lambda*x)))^a)^(b-1)
#' 
#' for x>0.
#' 
#' @return 
#' \code{dKMW} gives the density, \code{pKMW} gives the distribution 
#' function, \code{qKMW} gives the quantile function, \code{rKMW}
#' generates random deviates and \code{hKMW} gives the hazard function.
#' 
#' @export
#' @examples 
#' ## The probability density function  
#' curve(dKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), from = 0, to = 3, ylim = c(0, 2), col = "red", las = 1, ylab = "The probability density function")
#' 
#' ## The cumulative distribution and the Reliability function
#' par(mfrow = c(1, 2))
#' curve(pKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), from = 0, to = 3, col = "red", las = 1, ylab = "The cumulative distribution function")
#' curve(pKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2, lower.tail=FALSE), from = 0, to = 3, col = "red", las = 1, ylab = "The Reliability function")
#' 
#' ## The quantile function 
#' p <- seq(from = 0, to = 0.99999, length.out = 100)
#' plot(x = qKMW(p, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), y = p, xlab = "Quantile", las = 1, ylab = "Probability")
#' curve(pKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), from = 0, add = TRUE, col = "red")
#' 
#' ## The random function
#' hist(rKMW(n = 1000,alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), freq = FALSE, ylim = c(0, 2), xlab = "x", las = 1, main = "")
#' curve(dKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2),  from = 0, to = 3, add = TRUE, col = "red")
#' 
#' ## The Hazard function
#' curve(hKMW(x, alpha = 1, theta = 0.6, lambda = 2, a = 2, b = 1.2), from = 0, to = 1, ylim = c(0, 10), col = "red", ylab = "The Hazard function", las = 1)

dKMW<-function(x,alpha,theta,lambda,a,b, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  term <- alpha*(x^theta)*exp(lambda*x)
  loglik<- log(a*b*alpha) + (theta-1)*log(x) + log(theta+lambda*x) +
    (lambda*x-term) + (a-1)*log(1-exp(-term)) + 
    (b-1)*log(1-(1-exp(-term))^a)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}
#' @export
#' @rdname KMW
pKMW <- function(q,alpha,theta,lambda,a,b, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  cdf <- 1-(1-(1-exp(-alpha*(q^theta)*exp(lambda*q)))^a)^b 
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname KMW
qKMW <- function(p,alpha,theta,lambda,a,b, lower.tail = TRUE, log.p = FALSE){
  
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
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
  
  
  fda <- function(x,alpha,theta,lambda,a,b){
    1-(1-(1-exp(-alpha*(x^theta)*exp(lambda*x)))^a)^b 
  }
  
  fda1 <- function(x, alpha,theta,lambda,a,b, p) {fda(x, alpha,theta,lambda,a,b) - p}
  
  r_de_la_funcion <- function(alpha,theta,lambda,a,b, p) {
    uniroot(fda1, interval=c(0,1e+06), alpha,theta,lambda,a,b, p)$root
  }
  
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(alpha,theta,lambda,a,b, p)
  q
}
#' @export
#' @rdname KMW
rKMW <- function(n,alpha,theta,lambda,a,b){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qKMW(p,alpha,theta,lambda,a,b)
  r
}
#' @export
#' @rdname KMW
hKMW<-function(x,alpha,theta,lambda,a,b, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(alpha<=0 )) 
    stop(paste("alpha must be positive", "\n", ""))
  if (any(theta<=0)) 
    stop(paste("theta must be positive", "\n", ""))
  if (any(lambda<=0)) 
    stop(paste("lambda must be positive", "\n", ""))
  if (any(a<=0)) 
    stop(paste("a must be positive", "\n", ""))
  if (any(b<=0)) 
    stop(paste("b must be positive", "\n", ""))
  
  h <- dKMW(x,alpha,theta,lambda,a,b, log = FALSE)/pKMW(q=x,alpha,theta,lambda,a,b, lower.tail=FALSE, log.p = FALSE)
  h
}
