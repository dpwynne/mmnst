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

RDPCrossValidation = function (spikes, t.start = 0, t.end,
                               poss.lambda = seq(0, 10, by = 0.1),
                               max.J = 7, PSTH = FALSE,
                               max.diff = 0.005, pct.diff.plot = TRUE,
                               print.J.value = TRUE){
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

  T.data     <- t.end - t.start
  N.values   <- floor(2^(1:max.J))
  n.trials   <- length(spikes) #number of spike trains
  terminal.points <- lapply((t.end - t.start) / N.values, function(x) seq(t.start, t.end, x))
  ISE.matrix <- matrix(NA, nrow = length(poss.lambda), ncol = length(N.values))

  # Moved out of loop; was recomputing redundantly length(poss.lambda) times
  count.points = lapply(1:length(N.values), function(J) lapply(1:n.trials, function(i)
    # Use already built-in C function to do this
    #return(.Call(graphics:::C_BinCount, spikes[[i]], terminal.points[[J]], FALSE, FALSE)) # if this line replaces the next one, the code is two orders of magnitude faster, but it will mess up the R Check for CRAN because of the .Call function.
    #return(table(cut(spikes[[i]], terminal.points[[J]], right = FALSE, include.lowest = FALSE))) # this is slower version of the previous line
    return(Cpp_BinCount_sorted(spikes[[i]], terminal.points[[J]])) # this line is much faster the the previous one, bu still x2 slower than line 50
  ))

  if (PSTH){
    stop("Cross-validation based on PSTH is not yet implemented")  # Reza will fix this
  } else {

    for (J in 1:length(N.values)) {
      if (print.J.value)
        cat("J=", J, "\n")
      val <- N.values[J]
      for (lambda in poss.lambda) { #one pass = one entry in the lambda^th row of the ISE matrix
        ##for loop for f.hat.i
        f.hat.i <- matrix(0, nrow = n.trials, ncol = val)
        for (i in 1:n.trials) {
          xi <- spikes[[i]]
          if (length(xi) != 0)
            f.hat.i[i, ] = PoissonRDP(count.points[[J]][[i]], lambda) * (val / (T.data * length(xi)))
        }

        ##get f.hat estimate, square it, and integrate the squared estimate
        f.hat <- colMeans(f.hat.i)
        delta.t <- diff(terminal.points[[J]])
        integral.term <- sum(f.hat^2 * delta.t)

        # Some algebra to compute this matrix in one vectorized go
        f.hat.minus.i = matrix(f.hat * (n.trials / (n.trials - 1)),
                               nrow = n.trials, ncol = val, byrow = TRUE) - f.hat.i / (n.trials - 1)

        # Core of the change: vapply is faster than sapply when result structure known
        f.hat.minus.i.bar = vapply(1:n.trials, function(i) {
          if (length(spikes[[i]]) > 0) {
            til = spikes[[i]]
            points = terminal.points[[J]]

            # Do this in a vectorized way, with a C++ helper
            # Could probably still go faster here by assuming ordered til
            f.hat.minus.i.ti =
              sum(til == head(points, 1)) * f.hat.minus.i[i, 1] +
              sum(til == tail(points, 1)) * f.hat.minus.i[i, ncol(f.hat.minus.i)] +
              sum(f.hat.minus.i[i, fast_first_index_vec_ordered_x(til[til > min(points) & til < max(points)], points)])

            return(f.hat.minus.i.ti / length(spikes[[i]]))
          } else return(0)}, 0)

        C.V.ISE <- mean(f.hat.minus.i.bar) * 2
        ISE.matrix[which(poss.lambda == lambda), J] <- (integral.term - C.V.ISE)
      }
    }
  }

  ISE.for.J.equal.to.0 = -1/T.data
  min.ISE.for.each.J <- c(ISE.for.J.equal.to.0, apply(ISE.matrix,
                                                      2, min))
  min.ISE.overall <- min(min.ISE.for.each.J)
  min.threshold <- min.ISE.overall * (1 - max.diff)
  J.within.threshold <- which(min.ISE.for.each.J < min.threshold)
  min.J.allowable <- min(J.within.threshold) - 1
  if (min.J.allowable == 0) {
    cv.output.list <- list(J.ISE = 0, lambda.ISE = Inf, ISE = ISE.matrix)
  }
  else {
    lambda.min.ISE.index <- which.min(ISE.matrix[, min.J.allowable])
    cv.output.list <- list(J.ISE = min.J.allowable, lambda.ISE = poss.lambda[lambda.min.ISE.index],
                           ISE = ISE.matrix)
  }
  if (pct.diff.plot) {
    J.seq <- seq(0, max.J)
    pct.diff <- (min.ISE.overall - min.ISE.for.each.J)/min.ISE.overall *
      100
    plot(J.seq, pct.diff, pch = 16, ylab = "% diff",
         xlab = "J", main = "", ylim = c(0, max.diff * 200))
    abline(h = max.diff * 100, lty = "dotted")
  }
  # The leftmost point under the horizontal line corresponds to the chosen value of J
  return(cv.output.list)
}
