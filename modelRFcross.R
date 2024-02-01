library(plyr)
library(dplyr)
library(pROC)
library(randomForest)
library(caret)
library(e1071)


modelRFcross = function(data,
                    var.imp,
                    niter,
                    folds,
                    list.index.cross,
                    Nimportant = 'all'){
  
  modelRF0 = function(m,
                      data,
                      var.imp,
                      niter,
                      folds,
                      list.index.cross,
                      Nimportant = 'all'){
    
    if (!is.data.frame(data)) data = as.data.frame(data)
    
    data.train.cross = na.omit(data[list.index.cross$training[m,],])
    data.test.cross = na.omit(data[list.index.cross$validation[m,],])
    if (class(var.imp) == 'list') var.imp.m = var.imp[[m]]
    if (class(var.imp) == 'data.frame') var.imp.m = var.imp
    if (length(var.imp.m[,1]) < 2){
                 OOB = MCC = AUC = F1 = NA
      } else {
    
    if (Nimportant != 'all' & length(var.imp.m) > Nimportant)  var.imp = var.imp.m[1:Nimportant,1]
    if (Nimportant == 'all')  var.imp = var.imp.m[,1]
    
    wn = table(data.train.cross[,1])[[1]]/(table(data.train.cross[,1])[[1]]+table(data.train.cross[,1])[[2]])
    rf.model = randomForest(x = data.train.cross[,as.character(var.imp)] , 
                             y = as.factor(data.train.cross[,1]), 
                             cutoff = c(wn, 1-wn),
                             ntree = 1000)
    OOB = tail(c(as.data.frame(rf.model$err.rate)[,1]), n=1)
    rf.model.test = predict(rf.model, data.test.cross[,as.character(var.imp)])
    confusionMatrix = caret::confusionMatrix(rf.model.test, 
                                             as.factor(data.test.cross[,1]),
                                           dnn = c("Prediction", "Reference"))
    TP = confusionMatrix[[2]][1]
    FN = confusionMatrix[[2]][2]
    FP = confusionMatrix[[2]][3]
    TN = confusionMatrix[[2]][4]
    # print(confusionMatrix)
    MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    as.numeric.factor = function(x) {as.numeric(levels(x))[x]}
    AUC = pROC::auc(response = as.factor(data.test.cross[,1]), predictor = as.numeric.factor(rf.model.test))
    F1 = confusionMatrix[[4]][[7]]
    }
    result = list(OOB, MCC, AUC, F1)
    return(result)
  }
  
  N = folds*niter
  model.each.iter = lapply(1:N, function(m) modelRF0(m, 
                                                     data = data,
                                                     niter = niter,
                                                     var.imp = var.imp,
                                                     folds = folds,
                                                     list.index.cross = list.index.cross,
                                                     Nimportant = Nimportant))
  if (is.null(sapply(1:N, function(m) unlist(model.each.iter[[m]][[2]]))) == FALSE){
    
    OOB.mean = signif(mean(sapply(1:N, function(m) unlist(model.each.iter[[m]][[1]])),na.rm = TRUE), digits = 4)
    mcc.mean = signif(mean(sapply(1:N, function(m) unlist(model.each.iter[[m]][[2]])),na.rm = TRUE), digits = 4)
    auc.mean = signif(mean(sapply(1:N, function(m) unlist(model.each.iter[[m]][[3]])),na.rm = TRUE), digits = 4)
    f1.mean = signif(mean(sapply(1:N, function(m) unlist(model.each.iter[[m]][[4]])),na.rm = TRUE), digits = 4)
    
    
    OOB.sd = signif(sd(sapply(1:N, function(m) unlist(model.each.iter[[m]][[1]])),na.rm = TRUE)/sqrt(N), digits = 2)
    mcc.sd = signif(sd(sapply(1:N, function(m) unlist(model.each.iter[[m]][[2]])),na.rm = TRUE)/sqrt(N), digits = 2)
    auc.sd = signif(sd(sapply(1:N, function(m) unlist(model.each.iter[[m]][[3]])),na.rm = TRUE)/sqrt(N), digits = 2)
    f1.sd = signif(sd(sapply(1:N, function(m) unlist(model.each.iter[[m]][[4]])),na.rm = TRUE)/sqrt(N), digits = 2)
    
    OOB.list = sapply(1:N, function(m) unlist(model.each.iter[[m]][[1]]))
    mcc.list = sapply(1:N, function(m) unlist(model.each.iter[[m]][[2]]))
    auc.list = sapply(1:N, function(m) unlist(model.each.iter[[m]][[3]]))
    f1.list = sapply(1:N, function(m) unlist(model.each.iter[[m]][[4]]))
    
    
  } else {
    OOB.mean = mcc.mean = auc.mean = f1.mean = NA
    OOB.sd = mcc.sd = auc.sd = f1.sd = NA
    OOB.list = mcc.list = auc.list = f1.list = NA
  }
  
  result.model.frame =t(data.frame(OOB.mean, OOB.sd, 
                                    mcc.mean, mcc.sd,
                                    auc.mean, auc.sd,
                                    f1.mean, f1.sd))
  return(result.model.frame)                                    
}
