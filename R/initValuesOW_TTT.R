#' Initial values and search region for Odd Weibull distribution
#' 
#' @author Jaime Mosquera Gutiérrez \email{jmosquerag@unal.edu.co}
#' @family init
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
#'             environment from which \code{initValuesOW_TTT} is called.
#' @param local_reg a list of control parameters for LOESS. See 
#'                  \code{\link{loess.options}}.
#' @param interpolation a list of control parameters for interpolation function. See 
#'                  \code{\link{interp.options}}.  
#' @param ... further arguments passed to 
#'            \code{\link[EstimationTools]{TTTE_Analytical}}.                  
#'                  
#' @details 
#' This function performs a non-parametric estimation of the empirical total 
#' time on test (TTT) plot. Then, this estimated curve can be used so as to 
#' get suggestions about initial values and the search region for parameters 
#' based on hazard shape associated to the  shape of empirical TTT plot.
#' 
#' @example examples/examples_initValuesOW_TTT.R
#' @importFrom gamlss gamlss
#' @importFrom stats terms predict na.omit formula 
#' @importFrom survival is.Surv
#' @export                                                                                                                                               
initValuesOW_TTT <- function(formula, data=NULL,
                             local_reg = loess.options(),
                             interpolation = interp.options(), ...){
  if ( length(attr(terms(formula), "term.labels")) > 0 )
    warning('initValuesOW_TTT function only uses response variable.')
  mycall <- match.call()
  id_arg <- match(c('formula', 'data'), names(mycall),
                  nomatch=0)
  temp <- mycall[c(1L, id_arg)]
  temp[[1L]] <- quote(stats::model.frame)
  modfrm <- eval.parent(temp)
  y <- stats::model.extract(modfrm, 'response')
  outs <- EstimationTools::fo_and_data(y, formula, model_frame=modfrm, 
                                       data, fo2Surv = FALSE)
  fo <- outs$fo; data <- outs$data
  
  method <- if ( is.Surv(y) ){'censored'} else {'Barlow'}
  
  g1 <- EstimationTools::TTTE_Analytical(formula=fo, response=NULL,
                                         data=data, method=method,
                                         scale=TRUE, ...)
  g2 <- cbind(g1$`i/n`, g1$phi_n)
  g3 <- do.call("loess", list(g2[,2] ~ g2[,1], local_reg))
  g4 <- do.call(interpolation$interp.fun, list(x = g2[,1], y=predict(g3),
                                               interpolation$passing_args))
  
  dTTT_dp <- g4(seq(0,1,length.out = interpolation$length.out), deriv=1)
  d2TTT_dp2 <- g4(seq(0,1,length.out = interpolation$length.out), deriv=2)
  
  target <- diff(sign(d2TTT_dp2))
  inflex <- which( target != 0 )
  diff_val <- try(target[inflex], silent = TRUE)
  
  if ( length(inflex) < 2 ){
    if ( length(inflex) > 0 ){
      if (diff_val == 2){
        # Unimodal hazard
        sigma <- 0.6
        nu <- 7
        sigma.valid <- "all(sigma < 1)"
        nu.valid <- "TRUE"
        hazard_type <- "Unimodal"
      }
      if (diff_val == -2){
        # Bathtub hazard
        sigma <- 5
        nu <- 0.1
        sigma.valid <- "all(sigma > 1)"
        nu.valid <- "all(nu < 1) & all(nu > 0)"
        hazard_type <- "Bathtub"
      }
    } else {
      sign_search <- which(sign(d2TTT_dp2) > 0)
      if (is.na(sum(sign_search))){ # negative secod derivative
        # Decreasing hazard
        sigma <- 0.2
        nu <- 2
        sigma.valid <- "all(sigma < 1)"
        nu.valid <- "TRUE"
        hazard_type <- "Decreasing"
      } else { # positive second derivative
        # Increasing hazard
        sigma <- 2
        nu <- 6
        sigma.valid <- "all(sigma > 1)"
        nu.valid <- "all(nu > 0)"
        hazard_type <- "Increasing"
      }
    }
  } else {
    stop("Please, choose another initial values for parameters.
         Visit 'OW distribution' vignette to get further information.")
  }

  output <- list(formula=formula, response=y, 
                 sigma.start=sigma, nu.start=nu,
                 sigma.valid=sigma.valid, nu.valid=nu.valid,
                 local_reg=g3, interpolation=g4, TTTplot=g2,
                 hazard_type=hazard_type)
  class(output) <- "initValOW"
  return(output)
}

valid.region <- function(param, valid.values, formula, data){
  Error_valid <- "Please, define ther argument 'valid.values'
  in the right way. Visit 'OW distribution' vignette for further information."
  
  type <- class(valid.values)
  space <- paste0("init.OW$", param, ".valid")
  
  case_auto <- function(){
    if ( valid.values == "auto" ){
      init.OW <- initValuesOW_TTT(formula, data)
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
    if ( class(param_eval) == "try-error") stop(Error_valid)
    if ( class(param_eval) != "character") stop(Error_valid)
    return(param_eval)
  }
  
  res.param <- switch(type,
                      "character" = case_auto(),
                      "numeric" = stop(Error_valid),
                      "list" = case_manual())
  
  fun.param <- paste0("function(", param, ") ", res.param)
  return(eval(parse(text = fun.param)))
}