mult.add.model<-function(Time,ct,f,w0,eta,gamma,J=log(length(ct),2),pic=TRUE){
terminal.points<-identify.terminal.points(min(Time),max(Time),J)
multiplicative<-sapply(Time,theta.m,f=f,w0=w0,eta=eta,gamma=gamma,terminal.points=terminal.points,ct=ct)
additive<-sapply(Time,theta.a,f=f,w0=w0,eta=eta,gamma=gamma,terminal.points=terminal.points,ct=ct)

if(pic){
par(mfrow=c(3,1))
##plot ct, multiplicative, additive on same axes
plot(ct,type='s')
plot(multiplicative~Time,type='l')
plot(additive~Time,type='l')
}
return(list(mult=multiplicative,add=additive))
}
