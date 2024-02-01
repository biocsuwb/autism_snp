library(plyr)
library(dplyr)
library(Boruta)

# Input arguments 'varImpBoruta()':
# data - data set, the first column is descriptor values
# var.imp - vector of relevant variables, the same variable in each iteration
# folds - fold of crossvalidation
# niter - number of repeats
# list.index.cross - boolean mask train and test set

varimpBoruta = function(data,niter,folds,list.index.cross,maxRuns0){
  

  varBoruta = function(m,niter,folds,list.index.cross,maxRuns0){
    
   if (!is.data.frame(data)) data = as.data.frame(data)
    data.train.cross = na.omit(data[list.index.cross$training[m,],])
    #wn = table(data.train.cross[,1])[[1]]/table(data.train.cross[,1])[[2]]
    result.boruta <- Boruta(x = data.train.cross[,-1], y = as.factor(data.train.cross[,1]), maxRuns = maxRuns0)
    var.imp = getSelectedAttributes(result.boruta, withTentative = TRUE)
    #final.boruta <- data.frame(result.boruta$finalDecision)
    return(var.imp)
    } 
    
  N = niter*folds
  var.each.iter = lapply(1:N, function(m) varBoruta(m, 
                                                  niter = niter, 
                                                  folds = folds, 
                                                  list.index.cross = list.index.cross,
                                                  maxRuns = maxRuns0))
  return(var.each.iter)
}
