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
  cat("--------------------------------------------------------------------\n")
  cat("Initial Values\n")
  cat(paste0("sigma = ", object$sigma.start, "\n"))
  cat(paste0("nu = ", object$nu.start, "\n"))
  cat("--------------------------------------------------------------------\n")
  cat("Search Regions\n")
  cat("For sigma: ")
  cat(paste0(object$sigma.valid, "\n"))
  cat("For nu: ")
  cat(paste0(object$nu.valid, "\n"))
  cat("--------------------------------------------------------------------\n")
  cat("Hazard shape: ")
  cat(object$hazard_type)
  if ( !is.null(object$warning) ) warning(object$warning)
}
