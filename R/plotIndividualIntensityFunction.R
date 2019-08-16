plot.individual.intensity.functions<-function(spikes,train.numbers,f.hat.list,w0.hat.list,K.list,best.models,model.criteria,Time,terminal.points,ct,intensity.function.length=(1+max(Time)*10^(3-ceiling(log10(diff(range(Time))))))){
#produces intensity estimate and plots for specific spike trains identified by train.numbers
##f.hat.list, w0.hat.list are lists of f.hat and w0.hat for each model
##K.list is the list of matrices of eta.hat,gama.hat,fit for each model
##best.models is the output from step 8
##model.criteria should be "AICc" or "BIC" - determine which model to plot, default is best model by BIC
##intensity.function.length is the number of points in the discretized intensity func

cat("Plotting Individual Intensity Function for Best Models\n")

library(ggplot2)

#selected.models<-unique(best.models)
#selected.names<-unique(model.names)
##these two should have the same length

##change to this
selected.models<-best.models$model.numbers
selected.names<-best.models$model.names

nmodels<-length(selected.models)

etas<-lapply(K.list,"[[",1)
gamas<-lapply(K.list,"[[",2)

int.estimate<-vector("list",length=nmodels)
names(int.estimate)<-names(selected.names)

#par(mfrow=c(nmodels,1))

#for (i in 1:nmodels){

if (substr(model.criteria,1,3)=="AIC"){
i <- 2
}else{ ##default is BIC, otherwise AICc
i <- 3
}

##this ensures that only AICc (i=2) and BIC (i=3) run

list.number<-selected.models[i]

theta.t<-matrix(NA,intensity.function.length,length(train.numbers))

Time.vector<-seq(min(Time),max(Time),length=intensity.function.length)

for(ntrains in 1:length(train.numbers)){

f.hat<-f.hat.list[[list.number]]
w0.hat<-w0.hat.list[[list.number]][train.numbers[ntrains],]
eta.hat<-etas[[list.number]][train.numbers[ntrains],]
gamma.hat<-gamas[[list.number]][train.numbers[ntrains],]

if (list.number%%2==1){
theta.t[,ntrains]<-y<-sapply(Time.vector,theta.m,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)
}else{
theta.t[,ntrains]<-y<-sapply(Time.vector,theta.a,f=f.hat,w0=w0.hat,eta=eta.hat,gamma=gamma.hat,terminal.points=terminal.points,ct=ct)
}

}##end ntrains for loop

Time.fn<-rep(Time.vector,length(train.numbers))
t.fn<-c(theta.t)


trial.id <- paste0("Trial #",train.numbers)
trial.id.fn<-as.factor(matrix(trial.id,nrow=intensity.function.length,ncol=length(trial.id),byrow=T))

int.fns<-data.frame(Trial = trial.id.fn, Time = Time.fn, Intensity = t.fn)


plot.basics<-ggplot(int.fns,aes(Time,Intensity))+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(20.828,10.668,25.908,20.828),"mm"),axis.text=element_text(size=12))
plot.intfn<-plot.basics+geom_step(color="black",direction="hv",size=0.5)
plot.arrayed<-plot.intfn+facet_wrap(~Trial, labeller=label_value)

print(plot.arrayed)

return(int.fns)
}