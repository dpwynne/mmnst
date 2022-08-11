#' Recursive dyadic partitioning cross-validation
#'
#' Optimizes a penalized log-likelihood to find the optimal number of partitions for the recursive dyadic partitioning.
#'
#' @importFrom graphics plot abline
#'
#' @param spikes a list of spike trains.
#' @param t.start the starting time of the recording window; the default value is 0.
#' @param t.end the ending time of the recording window.
#' @param poss.lambda a numeric vector containing a grid of penalty values.
#' @param max.J the maximum resolution of the dyadic partitioning used the estimate the piecewise constant intensity function \eqn{c(t)}.
#' @param PSTH if TRUE, performs leave-one-train-out cross-validation for the c(t) estimate based on PSTH data.
#' If FALSE, performs leave-one-spike-out cross-validation for the c(t) estimate from each individual train.
#' @param max.diff the maximum allowance for the integrated squared error (ISE) of a smaller model to deviate from the overall minimum ISE.
#' @param pct.diff.plot a logical value indicating whether to produce a plot of the percentage difference (above minimum ISE) vs. J.
#' @param print.J.value a logical value indicating whether to print off the J value at each step of the cross-validation or not.
#'
#' @return A list of length 3 is returned returned.
#' The first item in the list is the optimal partition depth as computed by ISE (\eqn{\lambda}).
#' The second item in the list is the optimal penalty term as corresponding to that partition depth (J).
#' The third item in the list is a matrix containing the ISE values for all combinations of partition depth and penalty term.
#'
#' @export

RDPCrossValidation = function(spikes, t.start = 0, t.end,
                                     poss.lambda = seq(0, 10, by = 0.1), max.J = 7,
                                     PSTH = FALSE,
                                     max.diff = 0.005, pct.diff.plot = TRUE  , print.J.value = TRUE){
  # spikes = list of vectors; each vector represents the spike train for one of many repeated trials
  # Time = Time vector; can be just the start and end times of the recording; should be the same for all spike trains
  # poss.lambda = a grid of penalty values; by default, 0 to 10 in increments of 0.1
  # max.J = maximum value for J, where J is the resolution of the partitioning of the time vector;
  # by default, max.J = 7 for a 10-second train, which would correspond to a ~80 ms resolution
  # Note that setting max.J to very high values (above 7) can lead to computational issues and the resolution will be below the refractory period for relatively short spike trains
  # max.diff = the allowance for integrated squared error to deviate from the overall minimum integrated squared error, such that a model with lower J is preferred
  # as long as the integrated squared error is within max.diff relative change of overall minimum; by default, max.diff = 0.005 corresponding to a 0.5% allowance
  # when max.diff is large enough, the model with no partitioning (constant intensity) will be chosen; this case is equivalent to infinite penalty and will return J = 0
  # pct.diff.plot produces a plot of % difference vs. J

  T.data<-t.end-t.start

  ##get N.values
  N.values<-floor(2^(1:max.J))

  terminal.points<-vector("list",length=max.J)

  by.terminal<-(t.end-t.start)/N.values

  for (j in 1:length(N.values)){
    terminal.points[[j]] <- seq(t.start,t.end,by=by.terminal[j])
  }

  ##integrated Kullback-Leibler and integrated squared error
  ISE.matrix<-matrix(NA,nrow=length(poss.lambda),ncol=length(N.values))

  ##IKL = D_KL(f|f.hat) - using Wikipedia notation for our formula
  #IKL.matrix<-matrix(NA,nrow=length(poss.lambda),ncol=length(N.values))


  n.trials<-length(spikes) ##number of spike trains

  if (PSTH){
    stop("Cross-validation based on PSTH is not yet implemented")  # Reza will fix this
  } else {

  for(J in 1:length(N.values)){ ##one pass = one column of ISE matrix
    if(print.J.value) cat("J=",J,"\n")
    val <- N.values[J]
    for(lambda in poss.lambda){ #one pass = one entry in the lambda^th row of the ISE matrix
        ##for loop for f.hat.i
        f.hat.i<-matrix(NA,nrow=n.trials,ncol=val)
        for (i in 1:n.trials){
          xi <- spikes[[i]]
          count.points<-numeric(val)
          for (ii in 1:val){
            count.points[ii]<-length(xi[xi>=terminal.points[[J]][ii] & xi<terminal.points[[J]][ii+1]])
          }
          if(sum(xi) == 0){
            f.hat.i[i,] = rep(0 , length(count.points))
          }else{
            f.hat.i[i,]<-(PoissonRDP(count.points,lambda))*(val)/(T.data*length(xi))
          }
        }

        ##get f.hat estimate, square it, and integrate the squared estimate
        f.hat <- apply(f.hat.i,2,mean)
        delta.t<-diff(terminal.points[[J]])
        integral.term <- sum(f.hat^2*delta.t)

        ##for loop for f.hat.minus.i
        f.hat.minus.i <- matrix(NA,nrow=n.trials,ncol=val)
        f.hat.minus.i.bar<-numeric(n.trials)*NA
        f.hat.minus.i.log.bar<-numeric(n.trials)*NA
        for (i in 1:n.trials){
          temp.f<-f.hat.i[-i,,drop=F]
          f.hat.minus.i[i,]<-apply(temp.f,2,mean)
          if (length(spikes[[i]])>0){
            f.hat.minus.i.ti<-sapply(spikes[[i]], LeaveOneOutDensityEstimate,
                                     points=terminal.points[[J]],i=i,f.hat.minus.i=f.hat.minus.i)
          }else{
            f.hat.minus.i.ti = rep(0 , dim(temp.f)[2])
          }

          #f.hat.minus.i.bar for use in ISE
          f.hat.minus.i.bar[i]<-mean(f.hat.minus.i.ti)

          #f.hat.minus.i.log.bar for use in IKL
          #f.hat.minus.i.IKL<-f.hat.minus.i.ti
          #f.hat.minus.i.IKL[which(f.hat.minus.i.IKL <= 1e-20)] <- 1e-20
          #f.hat.minus.i.log.bar[i]<-mean(log(f.hat.minus.i.IKL))
        }

        C.V.ISE <- mean(f.hat.minus.i.bar)*2

        #C.V.IKL <- -mean(f.hat.minus.i.log.bar)

        ##ISE
        ISE.matrix[which(poss.lambda==lambda),J]<-(integral.term-C.V.ISE)

        ##IKL
        #IKL.matrix[which(poss.lambda==lambda),J]<- C.V.IKL

      }
    }

  }


  ISE.for.J.equal.to.0 = -1/T.data
  #IKL.for.J.equal.to.0 = -log(1/T.data)

  min.ISE.for.each.J <- c(ISE.for.J.equal.to.0, apply(ISE.matrix, 2, min))
  min.ISE.overall <- min(min.ISE.for.each.J)

  # the code below identifies the minimum J - the original code, with a plot of % diff vs. J, is in the archive
  min.threshold <- min.ISE.overall*(1-max.diff)

  J.within.threshold <- which(min.ISE.for.each.J < min.threshold)
  # note: this is a vector possibly starting at 1, which corresponds to J = 0
  # e.g. if J.within.threshold = c(2, 3, 7), the corresponding J are 1, 2, 6

  min.J.allowable <- min(J.within.threshold)-1

  if (min.J.allowable == 0){
    cv.output.list <- list(J.ISE = 0, lambda.ISE = Inf, ISE = ISE.matrix)
  } else {
    lambda.min.ISE.index <- which.min(ISE.matrix[,min.J.allowable])
    cv.output.list <- list(J.ISE = min.J.allowable, lambda.ISE = poss.lambda[lambda.min.ISE.index], ISE = ISE.matrix)
  }

  if(pct.diff.plot){
    J.seq <- seq(0, max.J)
    pct.diff <- (min.ISE.overall - min.ISE.for.each.J)/min.ISE.overall*100
    plot(J.seq, pct.diff, pch = 16, ylab = "% diff", xlab = "J", main = "", ylim = c(0,max.diff*200))
    abline(h = max.diff*100, lty = "dotted")
  }
  # The leftmost point under the horizontal line corresponds to the chosen value of J

  return(cv.output.list)
}
