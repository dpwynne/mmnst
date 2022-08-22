#' Select most common best model
#'
#' Identifies the model chosen as "best" for the most spike trains.
#'
#' @param selected.models a factor vector containing names of the "best" models as chosen by some IC.
#' @param model.names a factor vector containing the names of all models fit to the data.
#'
#' @return The name of the overall "best" model as chosen by plurality vote over all spike trains with a clear best model.
#'
#' @noRd

SelectMostCommonBestModel <- function(selected.models, model.names){
# selected.models is a factor vector of "best" models chosen by AIC or AICc or BIC
# should be one entry per spike train, but some can be NA
# model.names is the levels of the factors

model.frequencies <- table(selected.models)

if (dim(model.frequencies) == 0){
 return(model.names[1])
}

most.common.best.model <- model.names[which.max(model.frequencies)]
## by definition if two tie for the max, which.max selects the earlier index
## with the way we define our model list, this means that the simpler model will be chosen if there are ties
## if two models with same number of parameters tie, the multiplicative one will be chosen

return(most.common.best.model)

}
