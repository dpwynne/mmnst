#' Generate model names
#'
#' Creates a vector of names identifying the 2K+1 models.
#'
#' @param K.length the number of models to fit. If K is the maximum number of frequency components, then `K.length` = 2K+1.
#'
#' @return A character vector containing model names.
#' @export

GenerateModelNames <- function(K.length){
## K.length is the length of the K.list - i.e. the number of models
## this assumes a structure of muliplicative k, additive k, multiplicative k+1, additive k+1, etc.

model.list <- vector("character", K.length)

for (K in 1:(K.length-1)){
if (K%%2 == 1){
model.type = "Multiplicative"
}else{
model.type = "Additive"
}
model.number<-1+((K-1)%/%2)

model.list[K] <- paste(model.type, as.character(model.number))
}

model.list[K.length] <- "Nonperiodic"

names(model.list) <- model.list

return(model.list)
}
