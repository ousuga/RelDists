#' Initial values and search region for Odd Weibull distribution
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' @family initValOW
#' 
#' @description 
#' This function can be used so as to get suggestions about initial values
#' and the search region for parameter estimation in \code{OW} distribution.
#' 
#' @param formula an object of class \code{\link{formula}} with the response on
#'                the left of an operator \code{~}. The right side must be 
#'                \code{1}.
#' @param data an optional data frame containing the response variables. If 
#'             data is not specified, the variables are taken from the 
#'             environment from which \code{initValuesOW} is called.
#' @param local_reg a list of control parameters for LOESS. See 
#'                  \code{\link[EstimationTools]{loess.options}}.
#' @param interpolation a list of control parameters for interpolation function. See 
#'                  \code{\link[EstimationTools]{interp.options}}.  
#' @param ... further arguments passed to 
#'            \code{\link[EstimationTools]{TTTE_Analytical}}.                  
#'                  
#' @details 
#' This function performs a non-parametric estimation of the empirical total 
#' time on test (TTT) plot. Then, this estimated curve can be used so as to 
#' get suggestions about initial values and the search region for parameters 
#' based on hazard shape associated to the  shape of empirical TTT plot.
#' 
#' @return
#' Returns an object of class \code{c("initValOW", "HazardShape")} containing:
#' 
#' \itemize{
#' \item \code{sigma.start} value for \eqn{sigma} parameter of OW distribution.  
#' \item \code{nu.start} value for \eqn{nu} parameter of OW distribution.
#' \item \code{sigma.valid} search region for \eqn{sigma} parameter of OW distribution.  
#' \item \code{nu.valid} search region for \eqn{nu} parameter of OW distribution.
#' \item \code{TTTplot} Total Time on Test transform computed from the data.
#' \item \code{hazard_type} shape of the hazard function determined from the TTT
#' plot.
#' }
#' 
#' @example examples/examples_initValuesOW_TTT.R
#' @importFrom gamlss gamlss
#' @importFrom stats terms predict na.omit formula 
#' @importFrom survival is.Surv
#' @importFrom EstimationTools formula2Surv TTTE_Analytical TTT_hazard_shape fo_and_data
#' @importFrom EstimationTools loess.options interp.options
#' @importFrom BBmisc is.error
#' @export                                                                                                                                               
initValuesOW <- function(formula, data=NULL,
                             local_reg = loess.options(),
                             interpolation = interp.options(), ...){
  if ( length(attr(terms(formula), "term.labels")) > 0 )
    stop('initValuesOW function only uses response variable.')
  mycall <- match.call()
  id_arg <- match(c('formula', 'data'), names(mycall),
                  nomatch=0)
  temp <- mycall[c(1L, id_arg)]
  temp[[1L]] <- quote(stats::model.frame)
  modfrm <- eval.parent(temp)
  y <- stats::model.extract(modfrm, 'response')
  outs <- fo_and_data(y, formula, model_frame=modfrm,
                      data, fo2Surv=FALSE)
  fo <- outs$fo; data <- outs$data
  # if ( is.Surv(y) ){
  #   method <- 'censored'
  #   data <- as.data.frame(as.matrix(y))
  # } else {
  #   method <- 'Barlow'
  #   data <- as.data.frame(modfrm)
  # }
  
  dots <- substitute(...())
  args_matches <- match(names(formals(EstimationTools::TTT_hazard_shape)), 
                        names(dots), nomatch = 0)
  TTTE_params <- dots[args_matches]
  TTTE_dots <- dots[-args_matches]
  TTTE_dots <- if ( length(TTTE_dots) == 0 ){ NULL }
  
  Hazard_Shape <- do.call("TTT_hazard_shape",
                          args = c(list(formula = fo, data = data,
                                        local_reg = local_reg,
                                        interpolation = interpolation),
                                   TTTE_dots))
  formula <- Hazard_Shape$formula; y <- Hazard_Shape$response
  g3 <- Hazard_Shape$local_reg; g4 <- Hazard_Shape$interpolation
  g2 <- Hazard_Shape$TTTplot; hazard_type <- Hazard_Shape$hazard_type
  the_warning <- Hazard_Shape$warning

  if (!is.null(the_warning)){
    the_warning <- paste0("Non-parametric estimate for Empirical TTT",
                            " is irregular.\nPlease, ",
                            "use the 'plot()' method to see the TTT ",
                            "shape and set the search region manually in ",
                            "'gamlss()' if there is no conincidence between ",
                            "'Hazard_Shape()' and 'plot()'. Visit ",
                            "'OW distribution' vignette to get further ",
                            "information.")
  }
  
  if (is.error(g3) | is.nan(g3$s)){
    sigma <- NA;  nu <- NA; g3 <- NA
    sigma.valid <- NA; nu.valid <- NA
    warning(paste0("Problem with LOESS estimation. The sample",
                   "size may be too small"))
  } else {
    if ( !is.na(hazard_type) ){
      if (hazard_type == "Unimodal"){
        sigma <- 0.6
        nu <- 7
        sigma.valid <- "all(sigma < 1)"
        nu.valid <- "all(nu > 1/sigma)"
      }
      
      if (hazard_type == "Bathtub"){
        sigma <- 5
        nu <- 0.1
        sigma.valid <- "all(sigma > 1)"
        nu.valid <- "all(nu < 1/sigma)"
      }
      
      if (hazard_type == "Decreasing"){
        sigma <- 0.2
        nu <- 2
        sigma.valid <- "all(sigma < 1)"
        nu.valid <- "all(nu < 1/sigma)"
      }
      
      if (hazard_type == "Increasing"){
        sigma <- 2
        nu <- 6
        sigma.valid <- "all(sigma > 1)"
        nu.valid <- "all(nu > 1/sigma)"
      }
    } else {
      sigma <- NA;  nu <- NA;
      sigma.valid <- NA; nu.valid <- NA
    }
  }


  output <- list(formula=formula, response=y, strata=Hazard_Shape$strata, 
                 sigma.start=sigma, nu.start=nu,
                 sigma.valid=sigma.valid, nu.valid=nu.valid,
                 local_reg=g3, interpolation=g4, TTTplot=g2, 
                 hazard_type=hazard_type, warning=the_warning)
  class(output) <- c("initValOW", "HazardShape")
  return(output)
}