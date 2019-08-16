get.intensity.estimates<-function(spikes,Time,terminal.points,ct,freqrange,max.K){
##call this function after Time, terminal.points, and ct have been identified
##freqrange is a list of vectors of minimum/maximum frequencies to search over
##max.K is the number of periodic components in the most complex possible model

setup.pars<-setup.likelihoods(terminal.points)

##find all "peak" frequencies in the spike train
f.common.table<-find.top.freqs(spikes,Time,freqrange)

##if we want to have more frequencies than appear in the table, this cuts it down
max.K<-min(max.K,length(f.common.table))

f.hat.list<-w0.hat.list<-K.list<-vector("list",length=2*max.K)

##fit all possible models automatically
for (K in 1:max.K){
##position in the list of multiplicative and additive models for that K
nmodel.mult<-2*(K-1)+1
nmodel.add<-2*K



f.hat.list[[nmodel.mult]]<-f.hat.list[[nmodel.add]]<-pick.top.K.freqs(f.common.table,K)
w0.hat.list[[nmodel.mult]]<-w0.hat.list[[nmodel.add]]<-phase.estimation(spikes,f.hat.list[[nmodel.mult]])
K.list[[nmodel.mult]]<-fit.mult.model(spikes,f.hat.list[[nmodel.mult]],w0.hat.list[[nmodel.mult]],setup.pars,terminal.points,ct)
K.list[[nmodel.add]]<-fit.add.model(spikes,f.hat.list[[nmodel.add]],w0.hat.list[[nmodel.add]],setup.pars,terminal.points,ct)
}

bestModels<-model.select(K.list)

##comment this out for now - ellipses may not be usable
#plot.eta.gamma.ellipses(K.list,bestModels,f.hat.list)

##Estimates 1 averages parameters, Estimates 2 averages function
dev.new()
model.estimates.1<-extract.pars(f.hat.list,w0.hat.list,K.list,bestModels$model.numbers,bestModels$model.names,Time,terminal.points,ct)
dev.new()
model.estimates.2<-plot.theta.t(spikes,f.hat.list,w0.hat.list,K.list,bestModels$model.numbers,bestModels$model.names,Time,terminal.points,ct)

bestModels$model.names<-gsub(" ",".",bestModels$model.names)

names(model.estimates.1$estimate)<-unique(bestModels$model.names)
names(model.estimates.1$params)<-unique(bestModels$model.names)
names(model.estimates.2)<-unique(bestModels$model.names)

return(list(estimates.1=model.estimates.1,estimates.2=model.estimates.2))
}