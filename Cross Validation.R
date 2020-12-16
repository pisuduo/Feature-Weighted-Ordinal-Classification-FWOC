#### Cross Validation functions for FWOC (Feature weighted ordinal classification) ####
#### Author: Ziyang Ma ####
#### updated on 12/16/2020 ####

## -------- Cross Validation function --------
cv_fwoc=function(X,y,nk,r_vec,lambda_vec,d,weight_meth="spearman",sel="none",classifier = "LDA"){
  #### Cross Validation for FWOC####
  ## Input: X: n x p, y: n-length vector, nk: number of folds
  ## d: number of dimensions needed, 
  ## weight_meth: rank correlation for calculate feature weights, default spearman's rank
  ## sel: selection criterion for tuning parameter
  ## none: original criterion; 1, 2, 3, options avaialable.
  ## calssifier: afterwards classifier for projected data, options: LDA, KNN
  ## Output: 
  ## cv.acc: a matrix with classification accuracy for each combination of tuning parameters
  ## cv.acc.sd: a matrix with std of classification accuracy for each combination of tuning parameters
  ## best.acc, best.acc.sd, best.acc.lambda, best.acc.r:
  ## --->highest mean cv acc with standard error, with selected best parameters
  ## cv.tau and correponding: similar as cv.acc, metirc is kendall's tau
  ## cv.wt_acc and correponding: similar as cv.acc, metirc is weighted accuracy
  ## cv.avg_recall and correponding: similar as cv.acc, metirc is average recall
  ## cv.avg_prec and correponding: similar as cv.acc, metirc is average precision
  ## cv.wcost and correponding: similar as cv.acc, metirc is weighted misclassification cost
  ## cv.nz and correponding: similar as cv.acc, metirc is number of non-zero features
  ## Output:
  ## A list of metrics

  J=length(unique(y)) ## number of classes
  
  ## Initialization
  ## create nk folds
  set.seed(123)
  folds=createFolds(factor(y), k = nk, list = FALSE)
  cv.acc=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.acc.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.acc)=lambda_vec; rownames(cv.acc)=r_vec
  
  cv.wt_acc=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.wt_acc.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.wt_acc)=lambda_vec; rownames(cv.wt_acc)=r_vec
  
  cv.avg_recall=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.avg_recall.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.avg_recall)=lambda_vec; rownames(cv.avg_recall)=r_vec
  
  cv.avg_prec=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.avg_prec.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.avg_prec)=lambda_vec; rownames(cv.avg_prec)=r_vec
  
  cv.tau=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.tau.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.tau)=lambda_vec; rownames(cv.tau)=r_vec
  
  cv.wcost=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.wcost.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.wcost)=lambda_vec; rownames(cv.wcost)=r_vec
  
  cv.nz=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  cv.nz.sd=matrix(999,nrow=length(r_vec),ncol=length(lambda_vec))
  colnames(cv.nz)=lambda_vec; rownames(cv.nz)=r_vec
  
  X_tr_ls=list(); Y_tr_ls=list(); y_tr_ls=list()
  X_te_ls=list(); Y_te_ls=list(); y_te_ls=list()
  W_ls=list()
  
  
  for (k in 1:nk){
    te_ind=which(folds==k)
    
    ### reserve one fold for testing
    X_te_fold=X[te_ind,]
    y_te_fold=y[te_ind]
    Y_te_fold=create_dummyY(y_te_fold)
    
    ### select rest as the training data
    X_tr_fold=X[-te_ind,]
    y_tr_fold=y[-te_ind]
    Y_tr_fold=create_dummyY(y_tr_fold)
    
    
    X_tr_ls[[k]]=X_tr_fold; Y_tr_ls[[k]]=Y_tr_fold; y_tr_ls[[k]]=y_tr_fold
    X_te_ls[[k]]=X_te_fold; Y_te_ls[[k]]=Y_te_fold; y_te_ls[[k]]=y_te_fold
    
    ### calculating weights
    weight=abs(apply(X_tr_fold,2,cor,y=y_tr_fold,method=weight_meth))
    W=diag(1-weight,ncol=dim(X_tr_fold)[2],nrow=dim(X_tr_fold)[2])
    W_ls[[k]]=W
  }
  
  for(i in 1: length(r_vec)){
    print(i)
    for (j in 1:length(lambda_vec)){
      print(j)
      te_acc=rep(0,nk)
      wt_acc=rep(0,nk)
      avg_recall=rep(0,nk)
      avg_prec=rep(0,nk)
      te_tau=rep(0,nk)
      wcost=rep(0,nk)
      nz.swrlda=rep(0,nk)
      for (l in 1:nk){
        print(l)
        res=fwoc(X_tr_ls[[l]],Y_tr_ls[[l]],W_ls[[l]],r=r_vec[i],lambda = lambda_vec[j],k=d)
        res_ortho=ortho(res) ## orthogonalize the vectors
        
        nz.swrlda[l]=length(unique(c(unlist(apply(res_ortho,2,function(x) which(x!=0))))))
        
        X.train.proj=X_tr_ls[[l]]%*%res_ortho ### projection
        X.test.proj=X_te_ls[[l]]%*%res_ortho
        
        if (classifier == "LDA"){
          lda.out=lda(X.train.proj,grouping=as.factor(y_tr_ls[[l]]))
          lda.pred=predict(lda.out,X.test.proj)$class
          CM = confusionMatrix(lda.pred,reference = as.factor(y_te_ls[[l]]))
          CM=as.matrix(CM$table)
          te_acc[l]=mean(lda.pred==y_te_ls[[l]])
          te_tau[l]=cor(as.numeric(lda.pred),y_te_ls[[l]],method="kendall")
          wcost[l]=weightcost(y_te_ls[[l]],lda.pred)
          
        }
        if (classifier =="KNN"){
          knn.out_te=knn(X.train.proj,X.test.proj,y_tr_ls[[l]],k=5)
          te_acc[l]=mean(knn.out_te==y_te_ls[[l]])
          CM = confusionMatrix(knn.out_te,reference = as.factor(y_te_ls[[l]]))
          CM=as.matrix(CM$table)
          te_tau[l]=cor(as.numeric(knn.out_te),y_te_ls[[l]],method="kendall")
          wcost[l]=weightcost(y_te_ls[[l]],knn.out_te)
          
        }
        
        avg_recall[l] = mean(diag(CM)/colSums(CM))
        
        avg_prec[l]=avg_precison(CM)
        wt_acc[l] = mean(colSums(CM*AW(J)))
      }
      cv.acc[i,j]=mean(te_acc)
      cv.acc.sd[i,j]=sd(te_acc)
      cv.wt_acc[i,j]=mean(wt_acc)
      cv.wt_acc.sd[i,j]=sd(wt_acc)
      cv.avg_recall[i,j]=mean(avg_recall)
      cv.avg_recall.sd[i,j]=sd(avg_recall)
      cv.avg_prec[i,j]=mean(avg_prec)
      cv.avg_prec.sd[i,j]=sd(avg_prec)
      cv.tau[i,j]=mean(te_tau)
      cv.tau.sd[i,j]=sd(te_tau)
      cv.wcost[i,j]=mean(wcost)
      cv.wcost.sd[i,j]=sd(wcost)
      cv.nz[i,j]=mean(nz.swrlda)
      cv.nz.sd[i,j]=sd(nz.swrlda)
    }
  }
  if (sel=="none"){
    best.acc.ind=which(cv.acc == max(cv.acc), arr.ind = TRUE)
    best.tau.ind=which(cv.tau == max(cv.tau,na.rm = T), arr.ind = TRUE)
    best.wcost.ind=which(cv.wcost == min(cv.wcost), arr.ind = TRUE)
    
    best.acc=min(cv.acc)
    best.tau=max(cv.tau,na.rm = T)
    
    best.acc.lambda=lambda_vec[best.acc.ind[1,2]]
    best.acc.r=r_vec[best.acc.ind[1,1]]
    
    best.tau.lambda=lambda_vec[best.tau.ind[1,2]]
    best.tau.r=r_vec[best.tau.ind[1,1]]
    
    best.wcost.lambda=lambda_vec[best.wcost.ind[1,2]]
    best.wcost.r=r_vec[best.wcost.ind[1,1]]
  }
  if (sel=="1"){
    selection.acc=sel1_para_swrlda(cv.acc,"max")
    selection.tau=sel1_para_swrlda(cv.tau,"max")
    selection.wcost=sel1_para_swrlda(cv.wcost,"min")
    best.acc.lambda=lambda_vec[selection.acc[2]]; best.acc.r=r_vec[selection.acc[1]]
    best.tau.lambda=lambda_vec[selection.tau[2]]; best.tau.r=r_vec[selection.tau[1]]
    best.wcost.lambda=lambda_vec[selection.wcost[2]];best.wcost.r=r_vec[selection.wcost[1]]
  }
  if (sel=="2"){
    selection.acc=sel2_para_swrlda(cv.acc,"max")
    selection.tau=sel2_para_swrlda(cv.tau,"max")
    selection.wcost=sel2_para_swrlda(cv.wcost,"min")
    best.acc.lambda=lambda_vec[selection.acc[2]]; best.acc.r=r_vec[selection.acc[1]]
    best.tau.lambda=lambda_vec[selection.tau[2]]; best.tau.r=r_vec[selection.tau[1]]
    best.wcost.lambda=lambda_vec[selection.wcost[2]];best.wcost.r=r_vec[selection.wcost[1]]
  }
  if (sel=="3"){
    selection.acc=sel3_para_swrlda(cv.acc,"max")
    selection.tau=sel3_para_swrlda(cv.tau,"max")
    selection.wcost=sel3_para_swrlda(cv.wcost,"min")
    best.acc.lambda=lambda_vec[selection.acc[2]]; best.acc.r=r_vec[selection.acc[1]]
    best.tau.lambda=lambda_vec[selection.tau[2]]; best.tau.r=r_vec[selection.tau[1]]
    best.wcost.lambda=lambda_vec[selection.wcost[2]];best.wcost.r=r_vec[selection.wcost[1]]
  }
  
  
  
  return(list(cv.acc=cv.acc,cv.tau=cv.tau,
              cv.wt_acc =cv.wt_acc, cv.avg_recall =cv.avg_recall,cv.avg_prec =cv.avg_prec,
              cv.acc.sd=cv.acc.sd,cv.tau.sd=cv.tau.sd,
              cv.wcost=cv.wcost,cv.wcost.sd=cv.wcost.sd,
              cv.nz=cv.nz,cv.nz.sd=cv.nz.sd,
              best.acc.lambda=best.acc.lambda,best.acc.r=best.acc.r,
              best.tau.lambda=best.tau.lambda,best.tau.r=best.tau.r,
              best.wcost.lambda=best.wcost.lambda,best.wcost.r=best.wcost.r))
  
}