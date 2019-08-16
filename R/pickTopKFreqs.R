pick.top.K.freqs<-function(f.sorted,K, user.select = FALSE){
##pick the K frequencies that appear the most often
##in the sorted f.common table
cat("Identifying",K,"Most Common Peak Frequencies\n")
f.hat <- head(f.sorted,K)
print(f.hat) ##prints the top frequencies and the number of trials they exist at
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

cat(K,"Most Common Peak Frequencies Identified\n")
return(top.freqs)
}