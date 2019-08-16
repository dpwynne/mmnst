extract.pars<-function(f.hat.list,w0.hat.list,K.list,best.models,model.names,Time,terminal.points,ct){
##extracts parameters given set of best models,calls function that plots intensity estimate
##returns a list of time-intensity matrices and a list of parameter vectors
##f.hat.list, w0.hat.list are lists of f.hat and w0.hat for each model
##K.list is the list of matrices of eta.hat,gama.hat,fit for each model

cat("Determining Intensity Function for Best Models\n")

selected.models<-unique(best.models)
selected.names<-unique(model.names)
##these two should have the same length
nmodels<-length(selected.models)

etas<-lapply(K.list,"[[",1)
gamas<-lapply(K.list,"[[",2)

int.estimate<-params<-vector("list",length=nmodels)

par(mfrow=c(nmodels,1))

for (i in 1:nmodels){

list.number<-selected.models[i]

f.hat<-f.hat.list[[list.number]]
##for now, average; later, may do fancier math to get better estimate
w0.hat<-apply(w0.hat.list[[list.number]],2,mean)
eta.hat<-apply(etas[[list.number]],2,mean)
gamma.hat<-apply(gamas[[list.number]],2,mean)

f.hat.names<-paste0("f.",as.character(1:length(f.hat)))
w0.hat.names<-paste0("w0.",as.character(1:length(w0.hat)))
eta.hat.names<-paste0("eta.",as.character(1:length(eta.hat)))
gamma.hat.names<-paste0("gamma.",as.character(1:length(gamma.hat)))

params[[i]]<-c(f.hat,w0.hat,eta.hat,gamma.hat)
names(params[[i]])<-c(f.hat.names,w0.hat.names,eta.hat.names,gamma.hat.names)

int.estimate[[i]]<-plot.estimate(f.hat,w0.hat,eta.hat,gamma.hat,Time,terminal.points,ct,selected.names[i])

}
cat("Intensity Function Plotted\n")
return(list(estimate=int.estimate,params=params))
}

