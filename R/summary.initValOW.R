#' Summary of \code{initValOW} objects
#' 
#' @description 
#' This \code{summary} method displays initial values and search regions
#' for \code{\link{OW}} family.
#' 
#' @aliases summary.initValOW
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' 
#' @param object an object of class \code{initVal}, generated with 
#'               \code{\link{initValuesOW_TTT}}.
#' @param ... extra arguments
#'              
#' @export   
summary.initValOW <- function(object, ...){
  for (i in 1:length(object$strata)){
    cat("\n")
    cat(paste(object$type, names(object$strata)[i],'\n'))
    cat("--------------------------------------------------------------------\n")
    cat("Initial Values\n")
    cat(paste0("sigma = ", ifelse(is.null(object$sigma.start[[i]]), NA, 
                                  object$sigma.start[[i]]), "\n"))
    cat(paste0("nu = ", ifelse(is.null(object$nu.start[[i]]), NA, 
                               object$nu.start[[i]]), "\n"))
    cat("--------------------------------------------------------------------\n")
    cat("Search Regions\n")
    cat("For sigma: ")
    cat(paste0(ifelse(is.null(object$sigma.valid[[i]]), NA,
                      object$sigma.valid[[i]]), "\n"))
    cat("For nu: ")
    cat(paste0(ifelse(is.null(object$nu.valid[[i]]), NA, 
                      object$nu.valid[[i]]), "\n"))
    cat("--------------------------------------------------------------------\n")
    cat("Hazard shape: ")
    cat(ifelse(is.null(object$hazard_type[[i]]), NA, object$hazard_type[[i]]))
    cat("\n")
    cat("====================================================================\n")
  }
  if ( !is.null(object$warning) ) warning(object$warning)
}
