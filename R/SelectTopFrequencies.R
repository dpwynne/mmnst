#' Select top frequencies
#'
#' Identifies the K frequencies that appear most often in a given list of spike trains.
#'
#' @importFrom utils head
#'
#' @param f.sorted a sorted table. The names of the table, in order, are the most common frequencies in the periodograms of the individual spike trains.
#' This is usually output by the [FindTopFrequencies()] function.
#' @param K a scalar determining the number of periodic components to include in the model.
#' @param user.select whether to allow manual (user) control over accepting the identified frequencies. By default this is FALSE to allow running in batch mode.
#'
#' @return A numeric vector containing the K frequencies among the peaks of the (smoothed) periodogram which appear the most in the spike trains.
#'
#' @export

SelectTopFrequencies <- function(f.sorted,K, user.select = FALSE){
##pick the K frequencies that appear the most often
##in the sorted f.common table
cat("Identifying",K,"Most Common Peak Frequencies\n")
f.hat <- head(f.sorted,K)

for(i in 1:K){
  cat(paste(names(f.hat)[i], "Hz: ", f.hat[i], "trains"), "\n")
}

top.freqs<-as.numeric(names(f.hat))

if (user.select){  # if user.select, we want user to manually check and accept
while(1){
acceptable<-readline("Accept Frequencies? Type Y or N: ")
acceptable<-substr(acceptable,1,1)
if(acceptable == "Y" | acceptable == "y"){
break()
}
if(acceptable == "N" | acceptable == "n"){
cat("Identify Most Common Peak Frequencies Manually\n")
return()
}
cat("Response Not Accepted\n")
}  # end while loop
}

#cat(K,"Most Common Peak Frequencies Identified\n")
return(top.freqs)
}
