fit.all.mult.add.models <- function(K, spikes, f.common.table, setup.pars, terminal.points, ct, user.select = FALSE) {
## This abstracts Step 6 and 7 of the original script to run for each neuron
## For some value K, it fits Multiplicative 1-K, Additive 1-K, and the no periodicity model (i.e. c(t))

 f.hat.list<-list(c())
 w0.hat.list<-list(c())
 K.list<-list(c())

for(k in 1:K){
   f.hat.list[[2*k]] <- f.hat.list[[2*(k-1)+1]] <- pick.top.K.freqs(f.common.table, k, user.select = user.select)
   w0.hat.list[[2*(k-1)+1]] <-  w0.hat.list[[2*k]] <- phase.estimation(spikes,f.hat.list[[2*k]]) 
   K.list[[2*(k-1)+1]] <- fit.mult.model(spikes, f.hat.list[[2*(k-1)+1]], w0.hat.list[[2*(k-1)+1]], setup.pars, terminal.points, ct) 
   K.list[[2*k]] <- fit.add.model(spikes, f.hat.list[[2*k]], w0.hat.list[[2*k]], setup.pars, terminal.points, ct)
   names(f.hat.list)[c(2*k-1, 2*k)] <- paste(c("Multiplicative", "Additive"), k)
}

   f.hat.list[[2*K+1]] <- 0
   w0.hat.list[[2*K+1]] <-  phase.estimation(spikes,f.hat.list[[2*K+1]]) 
   K.list[[2*K+1]] <- fit.mult.model(spikes, f.hat.list[[2*K+1]], w0.hat.list[[2*K+1]], setup.pars, terminal.points, ct) 

   names(f.hat.list)[2*K+1] <- "Nonperiodic"
   names(w0.hat.list) <- names(f.hat.list)
   names(K.list) <- names(f.hat.list)

return(list(f.hat.list = f.hat.list, w0.hat.list = w0.hat.list, K.list = K.list))
}