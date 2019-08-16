check.fit<-function(l,k,n){
model.AIC<-(-2*l+2*k)
model.AICc<-(model.AIC+(2*(k+1)*(k+2)/(n-k-2)))
model.BIC<-(-2*l+log(n)*k)
return(c(model.AIC,model.AICc,model.BIC,l))
}