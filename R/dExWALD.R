#' The Ex-Wald distribution
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Ex-Wald distribution
#' with parameter \eqn{\mu}, \eqn{\sigma} and \eqn{\nu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param nu vector of the nu parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to 
#' response time data: An example using functions for the S-PLUS package. 
#' Behavior Research Methods, Instruments, & Computers, 36, 678-694.
#'
#' @seealso \link{ExWALD}.
#'
#' @details
#' The Wald distribution with parameters \eqn{\mu}, \eqn{\sigma} 
#' and \eqn{\nu} has density given by
#'
#' \eqn{f(x |\mu, \sigma)= x+1}
#'
#' @return
#' \code{dExWALD} gives the density, \code{pExWALD} gives the distribution
#' function, \code{qExWALD} gives the quantile function, \code{rExWALD}
#' generates random deviates.
#'
#' @example  examples/examples_dExWALD.R
#'
#' @export
dExWALD <- function(x, mu=1.5, sigma=1.5, nu=2, log=FALSE) {
  if (any(mu<=0))
    stop ("mu must be positive")
  if (any(sigma<=0))
    stop ("sigma must be positive")
  if (any(nu<=0))
    stop ("nu must be positive")
  
  # Freddy's code
  k <- mu^2 - 2/nu
  if (k < 0) {
    part1 <- -log(nu) + mu*sigma - sigma^2/(2*x) - x*mu^2/2
    element1 <- sqrt(-x*k/2)
    element2 <- sigma/sqrt(2*x)
    element3 <- ifelse(element1 < 30 & element2 < 25, 
                       rew(element1, element2), 
                       0)
    element3 <- abs(element3)
    part2 <- log(element3)
    res <- part1 + part2
  }
  else {
    k <- sqrt(k)
    part1 <- -log(nu)
    part2 <- sigma*(mu-k)-x/nu
    part3 <- pWALD(x, mu=k, sigma=sigma, log.p=TRUE)
    res <- part1 + part2 + part3  # Expression (3)
  }
  
  if (log)
    return(res)
  else
    return(exp(res))
}
# Vectorizing
dExWALD <- Vectorize(dExWALD)
#' @export
#' @rdname dExWALD
pExWALD <- function(q, mu=1.5, sigma=1.5, nu=2, 
                    lower.tail = TRUE, log.p = FALSE) {
  if (any(mu<=0))
    stop ("mu must be positive")
  if (any(sigma<=0))
    stop ("sigma must be positive")
  if (any(nu<=0))
    stop ("nu must be positive")
  
  cdf <- pWALD(q, mu, sigma) - nu * dExWALD(q, mu, sigma, nu)
  if (lower.tail == TRUE) {
    cdf <- cdf
  }
  else {
    cdf <- 1 - cdf
  }
  if (log.p == FALSE){
    cdf <- cdf}
  else {cdf <- log(cdf)}
  return(cdf)
}
#' @export
#' @rdname dExWALD
#' @importFrom stats uniroot
qExWALD <- function(p, mu=1.5, sigma=1.5, nu=2) {
  if (any(mu<=0))
    stop ("mu must be positive")
  if (any(sigma<=0))
    stop ("sigma must be positive")
  if (any(nu<=0))
    stop ("nu must be positive")
  
  # Begin auxiliar function
  my_aux <- function(x, p, mu, sigma, nu) 
    pExWALD(x, mu, sigma, nu) - p
  # End auxiliar function
  
  uniroot(my_aux, c(0.001, 10000), tol=0.0001, 
          mu=mu, sigma=sigma, nu=nu, p=p)$root
}
# Vectorizing
qExWALD <- Vectorize(qExWALD)
#' @export
#' @rdname dExWALD
#' @importFrom stats rexp
rExWALD <- function(n, mu=1.5, sigma=1.5, nu=2) {
  if (any(mu<=0))
    stop ("mu must be positive")
  if (any(sigma<=0))
    stop ("sigma must be positive")
  if (any(nu<=0))
    stop ("nu must be positive")
  
  rWALD(n, mu, sigma) + rexp(n, 1/nu)
}
#' Auxiliar function for the Ex-Wald distribution
#' @description This function is an auxiliar function.
#' @param x a value for x.
#' @param y a value for y.
#' @param block a value.
#' @param tol is the tolerance.
#' @param maxseries maximum value for series.
#' @return returns a vector with starting values.
#' @keywords internal
#' @export
uandv <- function(x, y, firstblock=20, block=0,
                  tol=.Machine$double.eps^(2/3), maxseries=20) {
  twoxy <- 2*x*y; xsq <- x^2; iexpxsqpi <- 1/(pi*exp(xsq))
  sin2xy <- sin(twoxy); cos2xy <- cos(twoxy)
  nmat <- matrix(rep((1:firstblock), each = length(x)), nrow = length(x))
  nsqmat <- nmat^2; ny <- nmat*y; twoxcoshny <- 2*x*cosh(ny)
  nsinhny <- nmat*sinh(ny); nsqfrac <- (exp(-nsqmat/4)/(nsqmat+4*xsq))
  u <- (2*pnorm(x*sqrt(2))-1)+iexpxsqpi*(((1-cos2xy)/(2*x))+2*
                                           ((nsqfrac*(2*x-twoxcoshny*cos2xy+nsinhny*sin2xy))%*%rep(1, firstblock)))
  v <- iexpxsqpi*((sin2xy/(2*x))+2*((nsqfrac*(twoxcoshny*sin2xy+
                                               nsinhny*cos2xy))%*%rep(1,firstblock)))
  n <- firstblock; converged <- rep(F, length(x))
  repeat {
    if ((block<1)||(n>=maxseries)) break
    else {
      if ((n+block)>maxseries) block <- (maxseries-n)
      nmat <- matrix(rep((n+1):(n+block), each=sum(!converged)),
                    nrow=sum(!converged))
      nsq <- nmat^2; ny <- nmat*y[!converged];
      twoxcoshny <- 2*x[!converged]*cosh(ny); nsinhny <- nmat*sinh(ny)
      nsqfrac <- (exp(-nsq/4)/(nsq+4*xsq[!converged]))
      du <- iexpxsqpi[!converged]*((2*nsqfrac*(2*x[!converged]-
                                                 twoxcoshny*cos2xy[!converged]+nsinhny*sin2xy[!converged]))
                                   %*%rep(1,block))
      dv <- iexpxsqpi[!converged]*((2*nsqfrac*(twoxcoshny*sin2xy[!converged]+
                                                nsinhny*cos2xy[!converged]))%*%rep(1, block))
      u[!converged] <- u[!converged]+du;
      v[!converged] <- v[!converged]+dv
      converged[!converged] <- ((du<tol)&(dv<tol))
      if (all(converged)) break
    }
  }
  cbind(u,v)
}

