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
#' @param initVal An \code{initValOW} object generated with \code{\link{initValuesOW_TTT}}
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
  link_funcs <- grep("*\\.link", family_call)
  
  # dist_type <- body(eval(family))[[5]][[2]][[2]]
  # cens_id <- grep(".*\\scens.*", dist_type)
  cens_id <- grep("cens.*", family_call)
  
  family <- if ( is.null(OWcall) ){family} 
  else {
    if (length(cens_id) > 0){
      eval(as.call( c(OWcall[[1]], OWcall[[2]])) )
    } else {
      if (length(OWcall)>1)
        eval(OWcall[[1]])
      else
        eval(OWcall)
    }
  }
  
  original_body <- body(family)
  size <- length(original_body)
  nopar <- body(family)[[5]][[2]]$nopar
  
  new_body <- original_body
  new_body[(size+1):(size+2)] <- c("","")
  
  new_body[1:(nopar+1)] <- original_body[1:(nopar+1)]
  new_body[[(nopar+2)]] <- NULL
  
  new_body[[(nopar+2)]] <- substitute(sigma.space <- 
                                        valid.region("sigma", valid.values, 
                                                     initVal), 
                                      list(valid.values=valid.values,
                                           initVal=initVal))
  
  new_body[[(nopar+3)]] <- substitute(nu.space <- 
                                        valid.region("nu", valid.values, 
                                                     initVal), 
                                      list(valid.values=valid.values,
                                           initVal=initVal))
  
  new_body[(nopar+4):(size+2)] <- original_body[(nopar+2):size]
  new_body[[(nopar+4)]][[2]]$sigma.valid <- substitute(sigma.space)
  new_body[[(nopar+4)]][[2]]$nu.valid <- substitute(nu.space)
  formals(family)$valid.values <- valid.values
  if ( length(link_funcs) > 0 ){
    if (length(cens_id) > 0){
      listOW <- as.list(as.call(as.list(OWcall[2])[[1]]))
    } else {
      listOW <- as.list(OWcall)
    }
    index <- 2:(length(OWcall))
    formals(family)[index] <- sapply(index, function(x) listOW[x])
  }
  body(family) <- new_body
  return(family)
}
#==============================================================================
# Selection of parametric space -----------------------------------------------
#==============================================================================
valid.region <- function(param, valid.values, initVal){
  Error_valid <- "Please, define ther argument 'valid.values'
  in the right way. Visit 'OW distribution' vignette for further information."
  
  type <- class(valid.values)
  space <- paste0("initVal$", param, ".valid")
  
  case_auto <- function(){
    if ( valid.values == "auto" ){
      param_space <- eval(parse(text = space))
    } else {
      stop(Error_valid)
    }
    return(param_space)
  }
  
  case_manual <- function(){
    list_pos <- paste0("valid.values$", param)
    param_eval <- try(eval(parse(text = list_pos)),
                      silent = TRUE)
    if ( inherits(param_eval, "try-error") ) stop(Error_valid)
    if ( !inherits(param_eval, "character") ) stop(Error_valid)
    return(param_eval)
  }
  
  res.param <- switch(type,
                      "character" = case_auto(),
                      "numeric" = stop(Error_valid),
                      "list" = case_manual())
  
  fun.param <- paste0("function(", param, ") ", res.param)
  return(eval(parse(text = fun.param)))
}