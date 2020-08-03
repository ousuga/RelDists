#' Mice mortality data
#'
#' The ages at death in weeks for male mice exposed to 240r of gamma radiation.
#' 
#' @docType data
#'
#' @usage data(mice)
#' 
#' @keywords datasets
#' 
#' @format A vector with 208 data points.
#' @importFrom graphics hist lines
#' @examples
#' data(mice)
#' hist(mice, main="", xlab="Time (weeks)", freq=FALSE)
#' lines(density(mice), col="blue", lwd=2)
"mice"