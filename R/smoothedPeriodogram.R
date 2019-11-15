smoothed.periodogram<-function(data,f,T.data,m=5,coef.step=0.01){
  #T.data<-(time.end-time.start); #m=max(1,floor(.1*length(data)))
  #commented out the above line, let's see if this works
  # In Brillinger's book, w is the angular frequency, but we use f: regular frequency
  function(f){
    # The next two lines have been kept this way to remind the difference
    # between angular and regular frequencies from Birllinger's book. The can be replaced by the
    # single line s.T<-round(f*T.data)
    w = 2*pi*f
    s.T<-round(w*T.data/2/pi)

    if (round(w%%pi,5)!=0){
      step.w<- -m:m; w.averaged.over<-2*pi*(s.T+step.w*coef.step)/T.data
      smoothed.perio<-(1/(2*m+1))*sum(periodogram(data , w.averaged.over/2/pi , T.data))
    }

    if (round(w%%(2*pi),5)==0 | (round(w%%(2*pi),5)==round(pi,5) & T.data%%2==0)){
      step.w<- 1:m; w.averaged.over<-w+2*pi*step.w*coef.step/T.data
      smoothed.perio<-(1/m)*sum(periodogram(data , w.averaged.over/2/pi , T.data))
    }

    if (round(w%%(2*pi),5)==round(pi,5) & T.data%%2!=0){
      step.w<- 1:m; w.averaged.over<-w-pi/T.data+2*pi*step.w*coef.step/T.data
      smoothed.perio<-(1/m)*sum(periodogram(data , w.averaged.over/2/pi , T.data))
    }
    return (smoothed.perio)
  }
}
