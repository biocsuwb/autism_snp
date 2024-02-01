
# Input arguments 'indexCross()':
# y - descriptor column
# folds - the fold of crossvalidation
# iterations -  the number of repetitions
# stratified - True/False
#
# Input arguments 'indexCross()':
# $training - the boolean mask for train subsets 
# $validation - the boolean mask for test subsets 


indexcross=function(y,folds,iterations,stratified) {
  
    stopifnot(folds>=1, length(y)>=folds, iterations>=1)
    folds=floor(folds)
        
    n.all=length(y)
    index=1:n.all
    validation=matrix(FALSE,iterations*folds,n.all)   
                
    if (stratified) {
         classes=unique(y)
         classes=classes[order(classes)]
         class.mask=y==t(matrix(classes,length(classes),length(y))) 
         n.class=colSums(class.mask) 
         samples=list()
         mask.class=matrix(0,n.all,length(classes))
         }
                                              
     for (it in 1:iterations) {
             if (stratified) for (cl in 1:length(classes)) samples[[cl]]=sample.int(n.class[cl]) else samples=sample.int(n.all)  
                                                            
        for (fold in 1:folds) {
              if (stratified) {
                  for (cl in 1:length(classes)) {
                           validation[(it-1)*folds+fold,]=validation[(it-1)*folds+fold,]|(index%in%which(class.mask[,cl])[samples[[cl]][(1:n.class[cl])>(fold-1)*n.class[cl]/folds & (1:n.class[cl])<=fold*n.class[cl]/folds]]) 
                  }
              } else { 
                validation[(it-1)*folds+fold,]=index%in%samples[(1:n.all)>(fold-1)*n.all/folds & (1:n.all)<=fold*n.all/folds]
              }
            }
    } 
    
    if (folds==1) training=validation else training=!validation
    return(list(training=training,validation=validation)) 
  } 