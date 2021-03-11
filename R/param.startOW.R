#' Initial values extraction for Odd Weibull distribution
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' 
#' @description 
#' This function can be used to extract initial values found with empirical 
#' time on test transform (TTT) through  \code{\link{initValuesOW_TTT}} function. 
#' This is used for parameter estimation in \code{OW} distribution.
#' 
#' @param param a character used to specify the parameter required. It can take the
#'              values \code{"sigma"} or \code{"nu"}.
#' @param initValOW an \code{initValOW} object generated with \code{\link{initValuesOW_TTT}}
#'                  function.
#'
#' @details
#' This function just gets initial values computed with \code{\link{initValuesOW_TTT}} 
#' for \code{OW} family. It must be called in \code{sigma.start} and \code{nu.start} 
#' arguments from \code{\link[gamlss]{gamlss}} function. This function is useful only
#' if user want to set start values automatically with TTT plot.
#' See example for an illustration.
#'     
#' @example examples/examples_param.startOW.R      
#' @export  
param.startOW <- function(param, initValOW){
  space <- paste0("initValOW$", param, ".start")
  res.param <- eval(parse(text = space))
  return(res.param)
}