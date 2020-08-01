#' Plot of \code{initValOW} objects
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' 
#' @description
#' Draws the empirical total time on test (TTT) plot and its non-parametric
#' (LOESS) estimated curve useful for searching initial values for \code{\link{OW}} 
#' family.
#' 
#' @aliases plot.initValOW
#' 
#' @param x an object of class \code{initVal}, generated with 
#'          \code{\link{initValuesOW_TTT}}.
#' @param xlab,ylab titles for x and y axes, as in \code{\link{plot}}.
#' @param curve_options a list with further arguments useful for customization
#'                      of non-parametric estimate plot.
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim the y limits (x1, x2) of the plot.
#' @param col The colors for lines and points. Multiple colors can be specified. This is 
#'            the usual color argument of \code{\link[graphics]{plot.default}}.
#' @param lty a vector of line types, see \code{\link{par}} for further information.
#' @param lwd a vector of line widths, see \code{\link{par}} for further information.
#' @param main a main title for the plot,
#' @param legend_options a list with fur further arguments useful for customization
#'                      of the legend of the plot. 
#' @param ... further arguments passed empirical TTT plot.
#'                      
#' @details 
#' This plot complements the use of \code{\link{initValuesOW_TTT}}. It is always 
#' advisable to use this function in order to check the result of non-parametric estimate 
#' of TTT plot. See the first example in \strong{Examples} section for an illustration.
#' 
#' The possible arguments for \code{...} can be consulted in 
#' \code{\link[graphics]{plot.default}}  and \code{\link{par}}.
#' 
#' @importFrom graphics par
#' @export   
plot.initValOW <- function(x, xlab="i/n", ylab=expression(phi(u)), xlim=c(0,1),
                         ylim=c(0,1), col = 1, lty=NULL, lwd=NA, main="", 
                         curve_options=list(col=2, lwd=2, lty=1),
                         legend_options=list(pos='top'),...){
  object <- x
  rm(x)
  par(xpd = TRUE, mar = par()$mar + c(0,0,0,7), par())
  plot(object$TTTplot[,1], object$TTTplot[,2], xlab=xlab, ylab=ylab, xlim=xlim, 
       ylim=ylim, main=main, col=col, lty=lty, lwd=lwd, ...)
  
  plot_options <- substitute(...())
  do.call("curve", c(list(expr=substitute(object$interpolation(x)), add=TRUE),
                     curve_options))
  
  legend_text <- c("Empirical TTT", "Spline curve")
  
  possible_pos <- c("top", "center", "bottom")
  if (!legend_options$pos %in% possible_pos) 
    stop((c("Select positions from the following list: \n \n",
            "  --> ",paste0(possible_pos, collapse=", "))))
  if (legend_options$pos == "center") legend_options$pos <- ""
  
  x <- paste0(legend_options$pos, "right")
  legend_options <- within(legend_options, rm(pos))
  legend_arguments <- c("y", "inset", "legend", "xpd", "col", "lty", "lwd")
  match_legend <- match(legend_arguments, legend_options, nomatch=0)
  match_legend <- which(match_legend != 0)
  
  if (length(match_legend) > 0) 
    stop(paste0("Argument(s)", legend_arguments[match_legend], "cannot be
                  manipulated. They have default unchangeable values."))
  
  do.call("legend", c(list(x, inset=c(-0.4,0), legend=legend_text,
                           pch=c(1,NA),
                           col=c(col, curve_options$col),
                           lty=c(lty,curve_options$lty), 
                           lwd=c(lwd,curve_options$lwd), xpd=TRUE),
                      legend_options))
}
