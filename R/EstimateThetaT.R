#' Estimate the intensity function \eqn{\theta(t)}
#'
#' Estimates average intensity function and computes bootstrap confidence intervals for a number of specified models.
#'
#' @param spikes a list of spike trains.
#' @param f.hat.list a list containing estimated frequencies for each model.
#' @param w0.hat.list a list containing estimated phases for each model.
#' @param K.list a list of matrices containing estimated eta and gamma parameters for each model and the estimated goodness-of-fit criteria (AIC, etc.) for each model.
#' @param models.to.fit a list containing the names of the models to fit an intensity function for and their corresponding indices in K.list.
#' Typically, this list contains either all models or only the models chosen by specific goodness-of-fit criteria.
#' @param t.start the starting time of the recording window; the default value is 0.
#' @param t.end the ending time of the recording window.
#' @param terminal.points a numeric vector containing the time points at which the piecewise constant estimate c(t) changes.
#' @param ct a numeric vector containing the estimated piecewise constant intensity function.
#' @param intensity.function.length the number of points in the discretized intensity function.
#' The larger this value is, the better the resolution.
#' For spike trains of 10 seconds or less, the default value corresponds to 10 ms resolution.
#' For spike trains of 10-100 seconds, the default value corresponds to 100 ms resolution.
#' @return  A list of length 2.
#' The first item in the list is a list of matrices containing the intensity estimates (average and bootstrap CI) for each model.
#' The second item in the list is a list of matrices containing the intensity estimates (for each spike train) for each model.
#'
#' @export
EstimateThetaT<-function(spikes, f.hat.list, w0.hat.list, K.list, models.to.fit, t.start = 0, t.end, terminal.points, ct, intensity.function.length=(1+(t.end - t.start)*10^(3-ceiling(log10(t.end - t.start))))){

  if(length(spikes) == 0){
    stop("No spike trains provided; cannot estimate intensity functions.")
  }

cat("Determining Intensity Function for Models\n")

selected.models<-models.to.fit[[2]] #$model.numbers
selected.names<-models.to.fit[[1]] #$model.names

nmodels<-length(selected.models)

etas<-lapply(K.list,"[[",1)
gamas<-lapply(K.list,"[[",2)

int.estimate<-vector("list",length=nmodels)
names(int.estimate)<-names(selected.names)

#par(mfrow=c(nmodels,1))

fitted.models<-c()

theta.t.list <- vector("list", length = nmodels)
#names(theta.t.list) <- c("AIC","AICc","BIC")

for (i in 1:nmodels){
#for (i in c(2,3)){
##this ensures that only AICc (i=2) and BIC (i=3) run

list.number<-selected.models[i]

theta.t<-matrix(NA,intensity.function.length,length(spikes))

Time.vector<-seq(t.start,t.end,length=intensity.function.length)

if (is.null(dim(ct))){ # if this is a vector of average ct, not a matrix of individual ct
  ct <- matrix(rep(ct, length(spikes)), nrow = length(spikes), byrow = TRUE)
}

f.hat<-f.hat.list[[list.number]]

#w0.hat<-apply(w0.hat.list[[list.number]], 2, mean)
# forces all phase to be average estimated phase over all spike trains
# This is not good when w0 estimates are near the boundaries of the window to average over

for(ntrains in 1:length(spikes)){

w0.hat<-w0.hat.list[[list.number]][ntrains,]
eta.hat<-etas[[list.number]][ntrains,]
gamma.hat<-gamas[[list.number]][ntrains,]

if (list.number%%2==1){
theta.t[,ntrains]<-sapply(Time.vector,ThetaMultiplicative,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct[ntrains,])
}else{
theta.t[,ntrains]<-sapply(Time.vector,ThetaAdditive,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct[ntrains,])
}

# these two lines are to properly scale theta(t) - this was rushed and should be revisited
theta.t.integral <- mean(theta.t[,ntrains])*(t.end - t.start) # approximation to integral of theta(t) using mean value theorem
theta.t[,ntrains] <-theta.t[,ntrains]/theta.t.integral*length(spikes[[ntrains]]) # scale theta.t such that it integrates to the number of spikes in the train

}##end ntrains for loop

theta.t.list[[i]] <- theta.t

theta.avg<-apply(theta.t,1,mean)


if (!(i %in% fitted.models)){
cat(paste("Determining Bootstrap Confidence Interval for",selected.names[i],"Model\n"))
theta.bootstrap <- BootstrappedTheta(theta.t)
theta.matrix <- cbind(Time.vector,theta.bootstrap[,1],theta.avg,theta.bootstrap[,2])
colnames(theta.matrix) <- c("Time","Lower","Average","Upper")
for (j in which(selected.models==selected.models[i])){
int.estimate[[j]]<-theta.matrix
fitted.models<-c(fitted.models,j)
}##end for
}##end if
}##end i for loop

names(theta.t.list) <- selected.names

cat("Intensity Function Estimated\n")
#return(int.estimate)
#return(int.estimate[c(2,3)])
#return(list(AIC = int.estimate[[1]],
#	AICc = int.estimate[[2]],
#	BIC = int.estimate[[3]],
#	individual.thetas = theta.t.list))
return(list(int.estimate = int.estimate,
	individual.thetas = theta.t.list))
}
