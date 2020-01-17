# It is important that the starting point of the widow to be 0 and the end point to be T
# in short, the data points tj must be smaller than T. Otherwise shif the data so that tj<T
# This needs documentation

periodogram<-function(data, f , T.data , centering = FALSE){ ## f is the frequency vector where the periodogram is evaluated.
  if(!centering){
    f.times.t<-f%*%t(data)
  }else{
    centralized.data = data - mean(data)
    f.times.t<-f%*%t(centralized.data)
  }
  return( (apply(cos(2*pi*f.times.t),1, sum)^2+apply(sin(2*pi*f.times.t),1,sum)^2)/(2*pi*T.data) )
}
