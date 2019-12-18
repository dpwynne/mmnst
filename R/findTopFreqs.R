#' Find highest-amplitude frequencies
#'
#' Finds the highest-amplitude frequency components in each spike train
#'
#' @importFrom utils tail
#'
#' @param spikes a list of spike trains
#' @param t.start the starting time of the recording window; the default value is 0
#' @param t.end the ending time of the recording window
#' @param freqrange a list of (non-overlapping) frequency ranges. Each item in the list should be a numeric vector of lowest and highest frequencies in the range.
#' @param q the number of highest-amplitude frequency components to find in each train
#' @param default.grid.spacing the spacing to use in the frequency search. This can be a single number reflecting the same grid spacing over all frequency ranges or a vector of the same length as freqrange
#' @param periodogram.window.size the number of points on each side of a given frequency to use when smoothing the periodogram
#' @param default.coef.step the coef.step value to pass to the smoothed periodogram
#'
#' @return a sorted table. The names of the table, in order, are the most common high-amplitude frequencies in the periodograms of the individual spike trains
#'
#' @export

FindTopFrequencies<-function(spikes, t.start = 0, t.end,
                         freqrange=list(c(2,30)), q=5,
                         default.grid.spacing = 1,
                         periodogram.window.size = 25,
                         default.coef.step = 0.01){
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

  ntrials<-length(spikes)
  T.data<-t.end-t.start

  ##get appropriate range of frequencies to check
  f = c()

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

    f<-c(f,seq(minfreq,maxfreq,by=max(default.grid.spacing[nrhythms],1/T.data)))

  }

  f.max<-vector("list",length=ntrials)
  for (replication.counter in 1:ntrials){
    x <- spikes[[replication.counter]]
    x <- x[(x >= t.start & x < t.end)] ##spikes within the time window

    if (length(x)>q){ ##we should have at least q+1 spikes if we are to estimate q frequencies

      smoothed.periodogram.f<-SmoothedPeriodogram(x,f,T.data, m = periodogram.window.size, coef.step = default.coef.step)

      estimate<-sapply(f,smoothed.periodogram.f)

      ##pick the q best frequencies
      q.max <- tail(sort(estimate),q)

      max.indx <- which(estimate %in% q.max)
      max.indx <- max.indx[max.indx>0]  # this line fixes a bug; we have forgotten what the bug is
      f.max[[replication.counter]] <- f[max.indx]
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
