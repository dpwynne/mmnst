identify.terminal.points<-function(t.start,t.end,J){
##t.start and t.end are the starting and ending times of the data collection
##J is identified from N*
terminal.points<-seq(t.start,t.end,length=(2^J+1)) ##vector of terminal points
}
