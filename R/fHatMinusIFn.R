f.hat.minus.i.fn <- function(til,i,points,f.hat.minus.i){

if (til < min(points) | til > max(points)) return(0)
if (til==min(points)) return(f.hat.minus.i[i,1])
if (til==max(points)) return(f.hat.minus.i[i,dim(f.hat.minus.i)[2]])
f.minus.i.t.index<-min(which(points>til))-1
return(f.hat.minus.i[i,f.minus.i.t.index])
}
