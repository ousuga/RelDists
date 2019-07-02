#' Fitting Different Parametric  \code{gamlss.family} Distributions (\code{RelDists} inclusive).
#'
#' @description 
#' This family of functions are a clones of \code{\link[gamlss]{fitDist}}, 
#' \code{\link{fitDistPred}}, \code{\link{chooseDist}} and
#' \code{\link{chooseDistPred}} of \code{\link[gamlss]{gamlss}} package. Type in
#' console \code{?gamlss::fitDist} for further information about the usage of this functions.
#' 
#' @param y the data vector
#' @param k a GAMLSS fitted model
#' @param type the type of distribution to be tried see details
#' @param try.gamlss this applies to functions fitDist() and fitDistPred(). It allows if gamlssML() fail to fit the model to try gamlss instead. This will slow up things for big data.
#' @param extra whether extra distributions should be tried, which are not in the type list
#' @param data the data frame where y can be found, only for functions fitDist() and fitDistPred()
#' @param trace whether to print during fitting. 
#' @param ... for extra arguments to be passed to gamlssML() to gamlss()
#' 
#' @details
#' The original \code{\link[gamlss]{fitDist}} has seven different family types; \code{realAll},
#' \code{realline, realplus, real0to1, counts} and \code{binom}. This implementation has an
#' additional function:
#' 
#' \itemize{
#'   \item \code{reldists}: distributions defined in the positive real line, implemented in
#'   \code{RelDists} package: "AddW", "BW", "EW", "ExW", "FWE", "GammaW",  "GMW", 
#'   "IW", "KW", "LW", "MW", "PL", "RW", "SZMW", "WG", "WGEE", "WP".
#' }
#' 
#' @seealso \code{\link[gamlss]{fitDist}}
#' 
#' @export 
fitDistPlus <- function(y,
                    k = 2, # for the AIC
                    type, 
                    try.gamlss = FALSE,  # whether to try the gamlss() if gamlssML() fails
                    extra = NULL,  # for extra distributions to include 
                    data = NULL, trace = FALSE, ...)
{
  
  reldists <- NULL
  if ("reldists" %in% type){
    reldists <- RelDists_families()
    output <- gamlss::fitDist(y = y, k = k, 
                              type = NULL, 
                              try.gamlss = try.gamlss, extra = c(reldists, extra), 
                              data = data, trace = trace, ...)
  } else{
    output <- gamlss::fitDist(y = y, k = k, 
                              type = type, 
                              try.gamlss = try.gamlss, extra = extra, 
                              data = data, trace = trace, ...)
  }
  return(output)
}
#' @importFrom stringr str_extract
RelDists_families <- function(){
  files_on_Rdir <- list.files()
  funcs <- str_extract(files_on_Rdir, "([A-Z]|[a-z])+")
  first_letters <- str_extract(files_on_Rdir, "([A-Z]|[a-z])")
  id_densities <- which(first_letters == "d")
  id_families <- which(first_letters != "d")
  
  # densities <- substring(funcs[id_densities], first = 2, last = 1000000L)
  families <- funcs[id_families]
  families <- families[families!="zzz"]
  families <- families[families!="fitDistPlus"]
  
  ourfamilies <- unique(families)
  return(ourfamilies)
}
