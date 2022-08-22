#' Information criterion-based multiscale model selection
#'
#' Selects the best model given a vector of AIC/AICc/BIC values for each model under consideration.
#'
#' @param IC a numeric vector of AIC, AICc, or BIC values for each model.
#' @param tolerance a very small positive scalar accounting for machine precision and roundoff errors in the optimization.
#' The simplest model within tolerance of the minimum IC value is selected as the best model.
#'
#' @return The index (in the IC vector) of the selected model, unless two equally simple models could be selected, in which case NA is returned.
#'
#' @noRd

SelectBestModelByIC <- function(IC, tolerance = 1e-6){
## selects the best model given a vector IC of AIC/AICc/BIC values for each model
## we are assuming the traditional mult-1, add-1, mult-2, add-2, ..., nonperiodic structure of the list
## tolerance is the amount by which we can go over the "minimum" IC value to select a smaller model
## tolerance should be a very small positive number, to account for precision of the machine and roundoff errors in the optimization

if (tolerance < 0) return(NA)  # if you have negative tolerance then something has gone wrong, return NA

IC.threshold <- min(IC)+tolerance

possible.index <- which(IC <= IC.threshold)

if (length(possible.index) == 1){  # if we have only one index
  return(possible.index)  # return that index
} else {
  max.K <- (length(IC) - 1)/2
  model.K.values <- c(rep(1:max.K, each = 2), 0)[possible.index]
  if (any(model.K.values == 0)){  # if we are comparing a nonperiodic model to one or more periodic models
    return(length(IC))  # return the index of the nonperiodic model
  } else {  # otherwise we are comparing only between periodic models
	  if( diff(model.K.values)[1] > 0) {  # if the two lowest-dimensional models have different dimensions
	    	return(possible.index[1])  # pick the lower-dimensional model
  	} else {
    		return(NA)  # otherwise we cannot decide between Mult-k and Add-k, return NA
  	}
  }
}

}
