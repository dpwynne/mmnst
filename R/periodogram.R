periodogram<-function(data, w, T.data){ ## w is the frequency vector where the periodogram is evaluated.
  w.data<-w%*%t(data)
  return( (apply(cos(w.data),1, sum)^2+apply(sin(w.data),1,sum)^2)/2/pi/T.data )
}