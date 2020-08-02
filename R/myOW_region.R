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
myOW_region <- function(family=OW, valid.values="auto"){
  
  previous.call <- sys.calls()
  match_all.calls <- sapply(previous.call,
                            function (x) match.call(gamlss, x))
  gamlss.pos <- which(regexpr("^?[g]amlss\\(formula", match_all.calls) == 1)
  
  if ( length(gamlss.pos) == 0 ){
    # Increasing hazard as default
    # sigma.space <- eval(parse(text = "function(sigma) all(sigma > 1)"))
    # nu.space <- eval(parse(text = "function(nu) all(nu > 0)"))
    # response <- NULL
    # fo <- NULL
    stop("This function must be called in gamlss() context.")
  } else {
    gamlss.call <- previous.call[[gamlss.pos]]
    call.est <- as.list(match.call(gamlss, gamlss.call))
    fo <- call.est$formula
    y_name <- all.vars(fo)[1]
    if (is.null(call.est$data) ){
      response <- as.name(y_name)
    } else {
      response <- call.est$data
    }
  }
  
  myOW_regioncall <- match.call()
  OWcall <- myOW_regioncall$family
  family_call <- deparse(OWcall)
  link_funcs <- grep("[.link]", family_call)
  
  family <- if ( is.null(OWcall) ){family} else {eval(OWcall[[1]])}
  original_body <- body(family)
  size <- length(original_body)
  nopar <- body(family)[[5]][[2]]$nopar
  
  new_body <- original_body
  new_body[(size+1):(size+3)] <- c("","","")
  
  new_body[1:(nopar+1)] <- original_body[1:(nopar+1)]
  
  new_body[[(nopar+2)]] <- quote(
    {
      y_name <- all.vars(fo)[1]
      if (!is.data.frame(response) ){
        y <- response
        data_frame <- data.frame(y)
        names(data_frame) <- y_name
      } else {
        data_frame <- response
        y <- paste0("data_frame", "$", y_name)
      }
      modfrm <-  stats::model.frame(fo, data=data_frame)
      gamlss.data <- EstimationTools::fo_and_data(y, fo, modfrm, 
                                                  data=data_frame, 
                                                  fo2Surv=FALSE)$data 
    }
  )
  new_body[[(nopar+3)]] <- quote(sigma.space <- 
                                      valid.region("sigma", valid.values, 
                                                   formula = fo, 
                                                   data = gamlss.data))
  new_body[[(nopar+4)]] <- quote(nu.space <- 
                                        valid.region("nu", valid.values, 
                                                     formula = fo,
                                                     data = gamlss.data))
  
  new_body[(nopar+5):(size+3)] <- original_body[(nopar+2):size]
  new_body[[(nopar+5)]][[2]]$sigma.valid <- substitute(sigma.space)
  new_body[[(nopar+5)]][[2]]$nu.valid <- substitute(nu.space)
  formals(family)$valid.values <- valid.values
  formals(family)$response <- response
  formals(family)$fo <- fo
  if ( length(link_funcs) > 0 ){
    index <- 2:(length(OWcall))
    formals(family)[index] <- sapply(index, function(x) as.list(OWcall)[x])
  } 
  body(family) <- new_body
  return(family)
}
