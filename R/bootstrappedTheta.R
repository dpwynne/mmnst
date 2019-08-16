bootstrapped.theta<-function(theta.t,B=1000){
##theta.t is a matrix with ntrains columns and length(time) rows
##B is the number of times to do bootstrap samples

ntrains<-dim(theta.t)[2]
bootstrap.matrix<-matrix(0,dim(theta.t)[1],B)

for (i in 1:B){
	sampled.intfns<-sample(1:ntrains,ntrains,replace=TRUE)
	new.sample.matrix<-theta.t[,sampled.intfns]
	new.mean<-apply(new.sample.matrix,1,mean)
	bootstrap.matrix[,i]<-new.mean
}##end for loop

bootstrapped.limits<-apply(bootstrap.matrix,1,quantile,probs=c(0.025,0.975))
if (dim(bootstrapped.limits)[2]>=dim(bootstrapped.limits)[1]) bootstrapped.limits<-t(bootstrapped.limits)
return(bootstrapped.limits)
}
