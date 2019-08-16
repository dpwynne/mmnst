setup.likelihoods<-function(terminal.points){
      D.i.plus.one<-terminal.points[-c(1)]
      D.i<-rev(rev(terminal.points)[-c(1)])
      DeltaDi=D.i.plus.one-D.i
	T.data<-max(terminal.points)-min(terminal.points)
	J<-log(length(terminal.points)-1,2)
return(list(DeltaDi=DeltaDi,Di.1=D.i.plus.one,Di.0=D.i,T.data=T.data,J=J))
}