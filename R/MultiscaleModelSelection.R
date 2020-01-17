#' Model selection for fitted multiscale models
#'
#' Finds the simplest multiscale model within a given tolerance of the "optimal" model as computed by an information criterion.
#'
#' @param K.list a list containing eta/gamma estimates and model fit criteria for each fitted model.
#'
#' @return A list of length 3.
#' The first item in the list is a character vector containing the names of the best model(s).
#' The second item in the list is a numeric vector containing the indices of the best model(s) in the list of models.
#' The third item in the list is a data frame containing the best model for each spike train analyzed by each fit criterion.
#'
#' @export

MultiscaleModelSelection <- function(K.list){

cat("Selecting Best Models\n")

#model.list<-vector(mode="character",length=length(K.list))

model.filled<-which(!sapply(K.list,is.null))

K.length<-length(K.list)

model.list <- names(K.list)

##fit is every third
#etas<-lapply(K.list,"[[",1)
#gamas<-lapply(K.list,"[[",2)
fits<-lapply(K.list,"[[",3)

ntrials<-dim(fits[[1]])[1]

fit.matrix<-matrix(NA,ntrials,3) ##AIC, AICc, BIC

rownames(fit.matrix)<-paste("Spike Train",seq(1,ntrials))
colnames(fit.matrix)<-c("AIC","AICc","BIC")

cat("Best Model by \t AIC \t AICc \t BIC\n")

#counterAIC <- counterAICc <- counterBIC <- rep(0,length(K.list))

for(itr in 1:ntrials){
AICall <- lapply(fits,function(x,it) x[it,1],it=itr)
AICcall <- lapply(fits,function(x,it) x[it,2],it=itr)
BICall <- lapply(fits,function(x,it) x[it,3],it=itr)
ll<-lapply(fits,function(x,it) x[it,4],it=itr)

criteriaAIC<-unlist(AICall)
criteriaAICc<-unlist(AICcall)
criteriaBIC<-unlist(BICall)

## AT THIS POINT WE HAVE ALL THE AIC'S, AICC'S, AND BIC'S FOR A SINGLE SPIKE TRAIN
## WHAT WE WANT TO DO WITH IT:
## 1. Find the minimum AIC, minimum AICc, minimum BIC
## 2. Find all indices in the AIC list within tolerance of minimum AIC
## 3. If the result of step 2 is one number, store that number
## 4. If the result of step 2 is multiple numbers, store all numbers
## 5. Do steps 2-4 for AICc and BIC
## 6. Once we have results for each spike train, first pass: count all the one-number rows
## 7. Second pass: In the multiple-number rows, for each row, compare the frequency of all indexes in that row
## 8. If two or more indices tie for most common, split among those numbers
## 9. Now there should be a total of length(spikes) "votes". If two indices are tied, select the lower index.
## 10. This means that multiplicative models will always be chosen over additive models with the same number of parameters.


## REVISED 4:30 PM on Tuesday 8/13/19
## Do Steps 1-5
## 6. Once have results for each spike train, only include the one-number rows.
## 7. Choose the best model based on only those rows. When the MLE is at the boundary of the closed space, the optimality properties
## may not hold [CHECK TO MAKE SURE THIS IS THEORETICALLY JUSTIFIABLE]

## REVISED 9:30 AM ON Wednesday 8/14/19
## Following revisions ending 11:30 PM Tuesday 8/13/19
## 1. Find minimum AIC/etc.
## 2. Find all models within tolerance of minimum AIC/etc.

indxAIC <- SelectBestModelByIC(criteriaAIC)
indxAICc <- SelectBestModelByIC(criteriaAICc)
indxBIC <- SelectBestModelByIC(criteriaBIC)

#indxAIC<-which.min(criteriaAIC)
#indxAICc<-which.min(criteriaAICc)
#indxBIC<-which.min(criteriaBIC)

#counterAIC[indxAIC] <- counterAIC[indxAIC]+1
#counterAICc[indxAICc]<- counterAICc[indxAICc]+1
#counterBIC[indxBIC] <- counterBIC[indxBIC]+1

## we should have NA if multiple models are best


#cat("Spike Train ",as.character(itr),"\t",model.list[model.filled[indxAIC]],"\t",model.list[model.filled[indxAICc]],model.list[model.filled[indxBIC]],"\n")

if (all(!is.na(c(indxAIC, indxAICc, indxBIC)))){
  fit.matrix[itr,]<-model.list[model.filled[c(indxAIC,indxAICc,indxBIC)]]
} else{  # if we have NA we have to go one by one
  which.na <- which(is.na(c(indxAIC, indxAICc, indxBIC)))
  which.not.na <- seq(1, 3)[-which.na]
  fit.matrix[itr,which.na] <- NA
  fit.matrix[itr, which.not.na] <- model.list[model.filled[which.not.na]]
 }

cat("Spike Train ",as.character(itr),"\t",fit.matrix[itr,1],"\t",fit.matrix[itr,2],"\t",fit.matrix[itr,3],"\n")

}



#cat("Warning: For some spike trains, sum of converged eta-hat values is 1.\nThis results in mulitiplicative and additive models being equivalent.")
## this means, periodic functions in the additive model become the same as the periodic functions in the multiplicative model,
## up to a multiplicative factor of c(t), i.e., the gamma-hat values in the additive model are the gamma-hat values of the multiplicative model times c(t)

level.names <- model.list

fit.output<-data.frame(AIC=factor(fit.matrix[,1],levels=level.names),AICc=factor(fit.matrix[,2],levels=level.names),BIC=factor(fit.matrix[,3],levels=level.names))

bestAIC <- SelectMostCommonBestModel(fit.output$AIC, model.names = level.names)
bestAICc <- SelectMostCommonBestModel(fit.output$AICc, model.names = level.names)
bestBIC <- SelectMostCommonBestModel(fit.output$BIC, model.names = level.names)

best.all <- c(bestAIC, bestAICc, bestBIC)

names(best.all)<-c("AIC","AICc","BIC")

model.numbers <- sapply(best.all, function(x) which(model.list %in% x))

cat("Best Model by: \n")
print(best.all,quote=F)
return(list(model.names=best.all,model.numbers=model.numbers,individual.fits=fit.output))
}

