find.top.freqs<-function(spikes,Time,freqrange=list(c(2,30)),q=5, default.grid.spacing = 1, periodogram.window.size = 25, default.coef.step = 0.01){
## spikes, a list of spike train data
## Time can be a vector of time points, a vector of start/end points,
## or the duration of recording (time.start assumed to be 0)
## we want to find the q highest-amplitude components in each train
## default.grid.spacing is the spacing to use in the frequency search, unless the spike trains are too short in time
## originally this was set to 0.2, but we are worried that 0.2 is too narrow to be biologically justifiable
## periodogram.window.size is the number of points on each side of a given frequency over which the periodogram is smoothed
## default.coef.step is the coef.step value to pass to smoothed periodogram
## i.e. the range of frequencies to look at to smooth at a frequency f is f +/- (periodogram.window.size)*(default.coef.step)
## and the spacing of the individual f values is default.grid.spacing

#cat("Finding Most Common Peak Frequencies\n")

if (length(Time) > 1){
time.start<-min(Time)
time.end<-max(Time)
}else{
time.start<-0
time.end<-Time
}
ntrials<-length(spikes)
T.data<-time.end-time.start

##get appropriate range of frequencies to check
w<-c()


if (length(default.grid.spacing) != length(freqrange)){

# if default grid spacing is one number, repeat for every range in freqrange
if (length(default.grid.spacing) == 1){
 default.grid.spacing <- rep(default.grid.spacing, length(freqrange))
} else {
 stop("There must be either one overall grid resolution or one grid resolution per interval of the frequency range to search over.")
}
}

for(nrhythms in 1:length(freqrange)){
minfreq<-freqrange[[nrhythms]][1]
maxfreq<-freqrange[[nrhythms]][2]

w<-c(w,2*pi*seq(minfreq,maxfreq,by=max(default.grid.spacing[nrhythms],1/T.data)))
}

f.max<-vector("list",length=ntrials)
for (replication.counter in 1:ntrials){
x <- spikes[[replication.counter]]
x <- x[(x >= time.start & x < time.end)] ##spikes within the time window

if (length(x)>q){ ##we should have at least q+1 spikes if we are to estimate q frequencies

smoothed.periodogram.w<-smoothed.periodogram(x,w,T.data, m = periodogram.window.size, coef.step = default.coef.step)

estimate<-sapply(w,smoothed.periodogram.w)

##pick the q best frequencies
q.max <- tail(sort(estimate),q)

max.indx <- which(estimate %in% q.max)
max.indx <- max.indx[max.indx>0]  # this line fixes a bug; we have forgotten what the bug is
f.max.i <- w[max.indx]/(2*pi)
f.max[[replication.counter]] <- f.max.i
}else{
f.max[[replication.counter]]<-rep(0,q)
}
#print(replication.counter) ##this just checks to make sure the for loop runs all the way through
}

##get the full list of frequencies
f.common <- unlist(f.max)
##summarize nonzero frequencies
f.common <- table(f.common[f.common>0])
f.sorted <- rev(sort(f.common)) ##ties broken by highest frequency
##f.sorted <- sort(f.common, decreasing=T) ##ties broken by lowest frequency
#cat("Most Common Peak Frequencies Found\n")
f.sorted.printable<-f.sorted
names(f.sorted.printable)<-paste0(names(f.sorted.printable),"Hz")
print(f.sorted.printable)
return(f.sorted)
}