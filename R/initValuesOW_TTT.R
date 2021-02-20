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
#' @importFrom EstimationTools formula2Surv TTTE_Analytical
#' @importFrom BBmisc is.error
#' @export                                                                                                                                               
initValuesOW_TTT <- function(formula, data=NULL,
                             local_reg = NULL,
                             interpolation = interp.options(), ...){
  if ( length(attr(terms(formula), "term.labels")) > 1 )
    stop('initValuesOW_TTT function only uses response variable.')
  mycall <- match.call()
  id_arg <- match(c('formula', 'data'), names(mycall),
                  nomatch=0)
  temp <- mycall[c(1L, id_arg)]
  temp[[1L]] <- quote(stats::model.frame)
  modfrm <- eval.parent(temp)
  y <- stats::model.extract(modfrm, 'response')
  outs <- fo_and_data_RelDists(y, formula, model_frame=modfrm, 
                               data, fo2Surv=FALSE)
  fo <- outs$fo; data <- outs$data
  
  method <- if ( is.Surv(y) ){'censored'} else {'Barlow'}
  
  dots <- substitute(...())
  if ( !('partition_method' %in% names(dots)) ){
    partition_method <- list(method = 'quantile-based',  folds = 3)
    dots[[(length(dots)+1)]] <- list(partition_method)
    names(dots[[length(dots)]])  <- 'partition_method'
  }
  
  args_matches <- match(names(formals(EstimationTools::TTTE_Analytical)), 
                        names(dots), nomatch = 0)
  TTTE_params <- dots[args_matches]
  TTTE_dots <- dots[-args_matches]
  TTTE_dots <- if ( length(TTTE_dots) == 0 ){ NULL }
  
  g1 <- do.call("TTTE_Analytical", 
                args = c(list(formula = fo, response = NULL,
                              data = data, method = method,
                              scale = TRUE), TTTE_params, TTTE_dots))
  
  TTTplot <- g3_list <- g4_list <- sigma <- nu <- sigma.valid <- 
    nu.valid <- hazard_type <- vector(mode = "list", 
                                      length = length(g1$strata))
  
  if ( is.null(local_reg) ) local_reg <- lapply(1:length(g1$strata), 
                                                function(i) loess.options())

  if ( length(local_reg) == 1 ) local_reg <- lapply(1:length(g1$strata),
                                                function(i) local_reg)
  the_warning <- NULL
  if ( all(grepl('(\\[.*\\)|\\[.*\\])', names(g1$strata))) ){
    type <- 'Interval'
  } else { type <- 'Level' }
  
  for (i in 1:length(g1$strata)){
    g2 <- cbind(g1$`i/n`[,i], g1$phi_n[,i])
    g2 <- na.omit(g2)
    TTTplot[[i]] <- g2
    
    g3 <- g3_list[[i]] <- try(do.call("loess", 
                                      c(list(formula=g2[,2] ~ g2[,1]),
                                        local_reg[[i]])), silent=TRUE)
    g4 <- g4_list[[i]] <- do.call(interpolation$interp.fun, 
                                  list(x = g2[,1], y=predict(g3),
                                       interpolation$passing_args))
    
    if (is.error(g3) | is.nan(g3$s)){
      sigma[[i]] <- NA;  nu[[i]] <- NA; g3 <- g3_list[[i]] <- NA
      sigma.valid[[i]] <- NA; nu.valid[[i]] <- NA
      hazard_type[[i]] <- NA
      warning(paste0("Problem with LOESS estimation. The sample",
                     "size may be too small"))
    } else {
      lout <- (length(y) - 1)*5
      dTTT_dp <- g4(seq(0,1,length.out = interpolation$length.out), deriv=1)
      d2TTT_dp2 <- g4(seq(0,1,length.out = interpolation$length.out), deriv=2)
      
      target <- diff(sign(d2TTT_dp2))
      inflex <- which( target != 0 )
      diff_val <- try(target[inflex], silent = TRUE)
      
      if ( length(inflex) < 2 ){
        if ( length(inflex) > 0 ){
          if (diff_val == 2){
            # Unimodal hazard
            sigma[[i]] <- 0.6
            nu[[i]] <- 7
            sigma.valid[[i]] <- "all(sigma < 1)"
            nu.valid[[i]] <- "all(nu > 1/sigma)"
            # nu.valid <- paste0("all(nu > 1/seq(1e-12, 1, length.out=",
            #                    as.character(as.name(lout)), "))")
            # all(nu > 1)"
            hazard_type[[i]] <- "Unimodal"
          }
          if (diff_val == -2){
            # Bathtub hazard
            sigma[[i]] <- 5
            nu[[i]] <- 0.1
            sigma.valid[[i]] <- "all(sigma > 1)"
            nu.valid[[i]] <- "all(nu < 1/sigma)"
            # nu.valid <- paste0("all(nu < 1/seq(1, 100, length.out=", 
            #                    as.character(as.name(lout)), "))")
            # "all(nu < 1) & all(nu > 0)"
            hazard_type[[i]] <- "Bathtub"
          }
        } else {
          sign_search <- any(sign(d2TTT_dp2) < 0) # if (is.na(sum(sign_search))){ 
          if (sign_search){# negative second derivative
            # Increasing hazard
            sigma[[i]] <- 2
            nu[[i]] <- 6
            sigma.valid[[i]] <- "all(sigma > 1)"
            nu.valid[[i]] <- "all(nu > 1/sigma)"
            # nu.valid <- paste0("all(nu > 1/seq(1,100, length.out=",
            #                    as.character(as.name(lout)), "))")
            # "all(nu > 0)"
            hazard_type[[i]] <- "Increasing"
          } else { # positive second derivative
            # Decreasing hazard
            sigma[[i]] <- 0.2
            nu[[i]] <- 2
            sigma.valid[[i]] <- "all(sigma < 1)"
            nu.valid[[i]] <- "all(nu < 1/sigma)"
            # nu.valid <- paste0("all(nu < 1/seq(1e-12, 1, length.out=",
            #                    as.character(as.name(lout)), "))")
            # "all(nu > 0)"
            hazard_type[[i]] <- "Decreasing"
          }
        }
      } else {
        the_warning <- paste0("Non-parametric estimate for ", 
                              paste(type, names(g1$strata)[i]),
                              " is irregular.\nPlease, ",
                              "use the 'plot()' method to see the TTT ", 
                              "shape and set the search region manually in ",
                              "'gamlss()' if there is no conincidence between ",
                              "'summary()' and 'plot()'. Visit ", 
                              "'OW distribution' vignette to get further ", 
                              "information.")
        warning(the_warning)
        criterion <- sapply(g2[,1], criteria, x_val=0, y_val=1, g3=g3)
        control1 <- all(criterion)
        # control2 <- isTRUE(criterion[1]) & isTRUE(criterion[length(g2[,1])])
        control2 <- all(criterion[2:(criterion[length(g2[,1])] - 1)])
        if ( control1 ){
          # Decreasing hazard
          sigma[[i]] <- 0.2
          nu[[i]] <- 2
          sigma.valid[[i]] <- "all(sigma < 1)"
          nu.valid[[i]] <- "all(nu < 1/sigma)"
          hazard_type[[i]] <- "Decreasing"
        } else if ( !control2 ){
          # Increasing hazard
          sigma[[i]] <- 2
          nu[[i]] <- 6
          sigma.valid[[i]] <- "all(sigma > 1)"
          nu.valid[[i]] <- "all(nu > 1/sigma)"
          hazard_type[[i]] <- "Increasing"
        } else {
          sigma[[i]] <- NA;  nu[[i]] <- NA;
          sigma.valid[[i]] <- NA; nu.valid[[i]] <- NA
          hazard_type[[i]] <- NA
        }
      }
    }
  }

  output <- list(formula=formula, response=y, predictor=g1$XB,
                 strata=g1$strata, sigma.start=sigma, nu.start=nu,
                 sigma.valid=sigma.valid, nu.valid=nu.valid,
                 local_reg=g3_list, interpolation=g4_list, type=type,
                 TTTplot=TTTplot, hazard_type=hazard_type, warning=the_warning)
  class(output) <- "initValOW"
  return(output)
}

valid.region <- function(param, valid.values, initVal, level){
  Error_valid <- "Please, define ther argument 'valid.values'
  in the right way. Visit 'OW distribution' vignette for further information."
  
  type <- class(valid.values)
  space <- paste0("initVal$", param, ".valid", "[[", level, "]]")
  
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
#==============================================================================
# Data preparation for TTT computation ----------------------------------------
#==============================================================================
#' @keywords internal
#'
fo_and_data_RelDists <- function(y, fo, model_frame, data, fo2Surv = TRUE){
  if ( !is.Surv(y) ){
    if ( fo2Surv ) fo <- EstimationTools::formula2Surv(model_frame)
    if ( missing(data) | is.null(data) ) data <- model_frame
  } else {
    if ( missing(data) | is.null(data) ){
      vars <- names(model_frame)
      ySurv <- vars[1L]
      yname <- gsub("Surv\\((.*?),.*", "\\1", ySurv)
      statusname <- gsub(paste0("Surv\\(", yname, ",(.*?)\\)"), "\\1", ySurv)
      right_hand <- attr(stats::terms(fo), 'term.labels')
      
      if (length(right_hand) == 0){
        factorname <- NULL
        data <- data.frame(y[,1], y[,2])
      } else {
        factorname <- as.character(right_hand[1])
        other_column <- model_frame[,2]
        data <- data.frame(y[,1], y[,2], other_column)
      }
      colnames(data) <- c(yname, statusname, factorname)
    }
  }
  return(list(data = data, fo = fo))
}
#==============================================================================
# Convexity criterion ---------------------------------------------------------
#
# f(lambda*x + (1 - lambda)*y) <= lambda*f(x) + (1 - lambda)*f(y)
#==============================================================================
#' @keywords internal
#'
criteria <- function(lambda, x_val, y_val, g3){
  f_xy <- predict(g3, newdata = c(x_val, y_val))
  right <- matrix(c(lambda, 1-lambda), ncol=2) %*% matrix(f_xy, nrow=2)
  right <- as.numeric(right)
  argument <- lambda*x_val + (1 - lambda)*y_val
  left <- predict(g3, newdata = argument)
  return(left <= right)
}