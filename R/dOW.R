#' @rdname OW
dOW<-function(x, mu, sigma, nu, log = FALSE){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  loglik<- log(mu) +log(sigma*nu) + (sigma-1)*log(x) +
    mu*(x^sigma) + (nu-1)*log(exp(mu*(x^sigma))-1) -
    2*log(1+(exp(mu*(x^sigma))-1)^nu)
  
  if (log == FALSE) 
    density<- exp(loglik)
  else 
    density <- loglik
  return(density)
}

#' @export
#' @rdname OW
pOW <- function(q,mu,sigma,nu, lower.tail=TRUE, log.p = FALSE){
  if (any(q<0)) 
    stop(paste("q must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  cdf <- 1 - (1 + (exp(mu*(q^sigma))-1)^nu )^(-1)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @export
#' @rdname OW
qOW <- function(p, mu, sigma, nu, lower.tail = TRUE, log.p = FALSE){
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- (1/mu)*(log( 1 + (-1+(1-p)^(-1))^(1/nu) ))^(1/sigma)
  q
}

#' @export
#' @rdname OW
rOW <- function(n, mu, sigma, nu){
  if(any(n<=0))
    stop(paste("n must be positive","\n",""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  n <- ceiling(n)
  p <- runif(n)
  r <- qOW(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname OW
hOW<-function(x,mu,sigma,nu){
  if (any(x<0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu<=0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma*nu<=0)) 
    stop(paste("Product sigma*nu must be positive", "\n", ""))
  
  h <- dOW(x, mu, sigma, nu, log = FALSE)/pOW(q = x, mu, sigma, nu, lower.tail=FALSE, log.p = FALSE)
  h
}