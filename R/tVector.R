tVector<-function(t.start,t.end,resolution){
##t.start and t.end are the starting and ending time points
##resolution is the resolution of the time signal
##either in time or frequency
if (resolution <= 0) return(c(t.start,t.end)) 
##non-positive resolution results in outputting start and end times only
if (resolution >= 1){
	tvec <- seq(t.start,t.end,by=(1/resolution))
}else{
	tvec <- seq(t.start,t.end,by=resolution)
}
return(tvec)
}
