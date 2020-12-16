#### Helper functions for FWOC (Feature weighted ordinal classification) ####
#### Author: Ziyang Ma ####
#### updated on 12/16/2020 ####

#---- Required packages ----

library(class)
library(caret) # split train and test -> createDataPartition
library(MASS)

#--- Center data matrix ----
center=function(X){
  ## Center the columns (variables) of a data matrix to zero mean
  ## Input: X: n x p data matrix
  n=dim(X)[1]
  m=colMeans(X)
  m.mat=matrix(rep(m,n),byrow=T,nrow=n)
  Xc=X-m.mat
  return(Xc)
}

#--- soft-thresholding function ----
Slambda=function(lambda,c){ 
  return(max(abs(c)-lambda,0)*sign(c))
}

#--- calculate feature weights (default: spearman) ----
calweights=function(X,y,meth='spearman'){
  ## calculate feature weights
  ## Inputs:
  ## X: n x p, data matrix
  ## y: n x 1, labels
  ## meth: method for calculating correlation, default: spearman, other option, kendall
  ## Outputs:
  ## w: p x 1, weight vector
  class_values=unique(y)
  J=length(class_values)
  p=dim(X)[2]
  
  if(meth=="discrete"){
    Cmean=matrix(0,nrow=J,ncol=p) ## calculate class means on each dimension
    for(i in 1:J){
      X_J=X[y==class_values[i],]
      Cmean[i,]=colMeans(X_J)
    }
    
    A_desc=matrix(0,nrow=J-1,ncol=J)
    
    for (i in 1:(J-1)){
      A_desc[i,i]=1
      A_desc[i,i+1]=-1
    }
    
    signsum=colSums(sign(A_desc%*%Cmean))
    w=ifelse(signsum==J-1|signsum==-(J-1),1,0)
    
  }else{
    w=matrix(0, ncol=1,nrow=p)
    for (i in 1:p){
      w[i]=abs(cor(X[,i],y,method=meth))
    }
    w=w/sum(w)
  }
  return (w)
} 

#--- create dummy Y matrix ----
create_dummyY=function(y){
  ## function to create dummy indicator Y matrix based on y
  ## Input: 
  ## y: n-length vector, labeled as 1 to K, K is the total number of classes
  n=length(y)
  K=length(unique(y))
  class_labels=sort(unique(y))
  Y=matrix(0,nrow=n,ncol=K)
  for (i in 1:n){
    for (j in 1:K){
      if (y[i]==class_labels[j]){
        Y[i,j]=1
      }
    }
  }
  return (Y)
}

#--- Calculating the between and within-covariance matrix ----
BW.mat1=function(X,Y){
  ## calcuating the between and within-covariance matrix
  ## Inputs:
  ## X: n x p design matrix
  ## Y: dummy indicator matrix
  ## Outputs:
  ## A list containing Sb (the between covariance matrix) and Sw (the within covariance matrix)
  Xc=center(X)
  
  n=dim(X)[1]
  J=dim(Y)[2]
  S11=t(Y)%*%Y/n
  S12=t(Y)%*%Xc/n
  S22=t(Xc)%*%Xc/n
  
  M=solve(S11)%*%S12
  Sb=t(M)%*%S11%*%M
  Sw=S22-Sb
  return(list(Sb=Sb,Sw=Sw))
  
}

# ---- Calculating the maximum lambda for tuning ----
l.max=function(X,Y,k){
  ## Input: X n x p
  ## Y: n x K dummy matrix
  ## k: number of disired discriminating vectors
  mats=BW.mat1(X,Y)
  Sw=mats$Sw
  Sb=mats$Sb
  A=Sb
  V=eigen(A)$vectors[,1:k]
  l.max=max(apply(V,1,function(x) sqrt(t(x)%*%x)))
  return(l.max)
}

#--- univariate ordinal logistic regression for feature selection----
uni.ordilg=function(X,y){
  ## univariate ordinal logistic regression for feature selection
  ## Input:
  ## X: n x p
  ## y: n-length vector
  ## Output: a vector of p-values for each feature, length of p
  
  y=as.factor(y)
  p=dim(X)[2]
  p.val=rep(0,p)
  names(p.val)=colnames(X)
  for (i in 1:p){
    print(i)
    mod=polr(y~X[,i],Hess=TRUE)
    ctable = coef(summary(mod))
    vals = pnorm(abs(ctable[,"t value"]),lower.tail = FALSE)*2
    p.val[i]=vals[1]
  }
  return(p.val)
}

####  ------------ Metrics ----------- ####
AW= function(J){
  ## calculate weight matrix to be muliplied with the confusion matrix ###
  ## Input: J: number of classes
  ## Output: WW: weight matrix
  WW = matrix(0,ncol=J,nrow=J)
  for(i in 1:J){
    for (j in 1:J){
      WW[i,j] =1-abs(i-j)/J
    }
  }
  WW = apply(WW,2,function(x){x/sum(x)})
  return(WW)
}

#--- weighted cost for misclassfication ----
weightcost=function(actual,pred,d=1){
  ### Input:
  ### actual: a vector of actual class memberships, numerical labels representing orders
  ### pred: a vector of predicted class memberships, numerical labels representing orders
  ### d: polynomial order of absolute cost, default=1
  cost=0
  actual=as.numeric(actual)
  pred=as.numeric(pred)
  if (length(actual)!= length(pred)){
    stop("Error: length of predictions do not match!")
  }
  else{
    cost=sum(abs(actual-pred)^d)/length(actual)
    return(cost)
  }
}

## ----average class precision----
avg_precison = function(CM){
  ## Input: CM, confusion matrix
  ## Output: average class precision
  m = dim(CM)[1]
  prec = rep(0,m)
  for (i in 1:m){
    if(rowSums(CM)[i]!=0){
      prec[i]=diag(CM)[i]/rowSums(CM)[i]
    }
  }
  return(mean(prec))
}

####----------- cross  validation tuning parameter selection criteria ------ ####

## the first selection criterion, in case there is a tie among tuning parameters
sel1_para_swrlda=function(metric,option="max"){
  ## metric:  a matrix with metric, each grid represents a combination of tuning parameters
  ## principle: 
  ## among all the "best" candidates, filter r not greater than median (r), 
  ## among these candidates, select lambda with largest lambda, 
  ## among the rest, select the minimum r.
  ## option="max" or "min"
  if (option=="max"){
    ind=which(metric==max(metric,na.rm=TRUE),arr.ind=TRUE)
  }else{
    ind=which(metric==min(metric,na.rm = TRUE),arr.ind=TRUE)
  }
  if (is.null(dim(ind))==TRUE){
    return(ind)
  }else{
    ## step 1, filter out r greater than median (r)
    ind1=ind[which(ind[,1]<=median(ind[,1])),]
    if (is.null(dim(ind1))==TRUE){
      return (ind1)
    }else{
      ## step 2, select largest lambda
      ind2=ind1[which(ind1[,2]==max(ind1[,2])),]
      if (is.null(dim(ind2))==TRUE){
        return(ind2)
      }else{
        ## step 3, among the rest, select minimum r
        final_ind=ind2[order(ind2[,1],decreasing=FALSE),][1,]
        return (final_ind)
      }
    }
  }
}
## the second selection criterion, smallest r then largest lambda
sel2_para_swrlda=function(metric,option="max"){
  ## metric:  a matrix with metric, each grid represents a combination of tuning parameters
  ## principle: 
  ## among all the "best" candidates, select smallest r
  ## among the rest, select largest lambda
  if (option=="max"){
    ind=which(metric==max(metric,na.rm=TRUE),arr.ind=TRUE)
  }else{
    ind=which(metric==min(metric,na.rm = TRUE),arr.ind=TRUE)
  }
  if (is.null(dim(ind))==TRUE){
    return(ind)
  }else{
    ## step 1, select smallest r
    ind1=ind[which(ind[,1]==min(ind[,1])),]
    if(is.null(dim(ind1))==TRUE){
      return (ind1)
    }else{
      ## step 2, select largest lambda
      ind2=ind1[which(ind1[,2]==max(ind1[,2])),]
      return (ind2)
    }
  }
}
## the third selection criterion, largest lambda then smallest r
sel3_para_swrlda=function(metric,option="max"){
  ## metric:  a matrix with metric, each grid represents a combination of tuning parameters
  ## principle: 
  ## among all the "best" candidates, select largest lambda
  ## among the rest, select  smallest r
  if (option=="max"){
    ind=which(metric==max(metric,na.rm=TRUE),arr.ind=TRUE)
  }else{
    ind=which(metric==min(metric,na.rm = TRUE),arr.ind=TRUE)
  }
  if (is.null(dim(ind))==TRUE){
    return(ind)
  }else{
    ## step 1, select largest lambda
    ind1=ind[which(ind[,2]==max(ind[,2])),]
    if (is.null(dim(ind1)[1])==TRUE){
      return (ind1)
    }else{
      ## step 2, select smallest r
      ind2=ind1[which(ind1[,1]==min(ind1[,1])),]
      return (ind2)
    }
  }
}

## the forth selection criterion, "Min-Min" principle
sel4_para_swrlda=function(metric,num_features){
  ## metric:  a matrix with metric, each grid represents a combination of tuning parameters
  ## num_features: a matrix with number of non-zero features, each grid represents a combination of tuning parameters
  ## principle: 
  ## among all the "best" candidates, select the pairs with the smallest non-zero features
  ## if replicates, select the first occur ones.
  if (option=="max"){
    ind=which(metric==max(metric,na.rm=TRUE),arr.ind=TRUE)
  }else{
    ind=which(metric==min(metric,na.rm = TRUE),arr.ind=TRUE)
  }
  if (is.null(dim(ind))==TRUE){
    return(ind)
  }else{
    ## step 1, select largest lambda
    num=num_features[ind]
    
    ind1=ind[which(ind[,2]==max(ind[,2])),]
    if (is.null(dim(ind1)[1])==TRUE){
      return (ind1)
    }else{
      ## step 2, select smallest r
      ind2=ind1[which(ind1[,1]==min(ind1[,1])),]
      return (ind2)
    }
  }
}


#### ----------- After Main function -----####
#--- find the orthonormal basis of matrix ----
ortho=function(mat){
  ## find the orthonormal basis of matrix
  ## Input: 
  ## mat: matrix
  ## Output:
  ## q: the orthonormal basis of mat
  qrres=qr(mat)
  q=qr.Q(qrres)
  return(q)
}
