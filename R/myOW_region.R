#' Custimized region search for odd Weibull distribution
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' 
#' @description 
#' This function can be used to modify \code{OW} \code{gamlss.family} object 
#' in order to  set a customized region search for \code{gamlss()} function.
#' 
#' @param family The \code{\link{OW}} family. This arguments allows the user to
#'               modify input arguments of the family, like the \code{link}
#'               functions.
#' @param valid.values a list of character elements specifying the region for 
#' \code{sigma} and/or \code{nu}. See \strong{Details} and \strong{Examples}
#' section to learn about its use.
#' @param initVal an \code{initValOW} object generated with \code{\link{initValuesOW_TTT}}
#' function.
#' 
#' @details 
#' This function was created to help users to fit \code{OW} distribution easily 
#' bounding the parametric space for \code{sigma} and \code{nu}.
#' 
#' The \code{valid.values} must be defined as a list of characters containing a call
#' of the \code{\link{all}} function.
#' 
#' @example examples/examples_myOW_region.R
#' @export
myOW_region <- function(family=OW, valid.values="auto", initVal){
  
  response <- initVal$response
  mydata <- paste0("data.frame(", all.vars(initVal$formula)[1], "=response)")
  
  myOW_regioncall <- match.call()
  OWcall <- myOW_regioncall$family
  family_call <- deparse(OWcall)
  link_funcs <- grep("[.link]", family_call)
  
  family <- if ( is.null(OWcall) ){family} else {eval(OWcall[[1]])}
  original_body <- body(family)
  size <- length(original_body)
  nopar <- body(family)[[5]][[2]]$nopar
  
  new_body <- original_body
  new_body[(size+1):(size+2)] <- c("","")
  
  new_body[1:(nopar+1)] <- original_body[1:(nopar+1)]

  new_body[[(nopar+2)]] <- substitute(sigma.space <- 
                                        valid.region("sigma", valid.values, 
                                                     formula, data), 
                                      list(valid.values=valid.values,
                                           formula=initVal$formula,
                                           data=eval(parse(text=mydata))))
  
  new_body[[(nopar+3)]] <- substitute(nu.space <- 
                                        valid.region("nu", valid.values, 
                                                     formula, data), 
                                      list(valid.values=valid.values,
                                           formula=initVal$formula,
                                           data=eval(parse(text=mydata))))
  
  new_body[(nopar+4):(size+2)] <- original_body[(nopar+2):size]
  new_body[[(nopar+4)]][[2]]$sigma.valid <- substitute(sigma.space)
  new_body[[(nopar+4)]][[2]]$nu.valid <- substitute(nu.space)
  formals(family)$valid.values <- valid.values
  if ( length(link_funcs) > 0 ){
    index <- 2:(length(OWcall))
    formals(family)[index] <- sapply(index, function(x) as.list(OWcall)[x])
  } 
  body(family) <- new_body
  return(family)
}
