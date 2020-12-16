#### Main functions for FWOC (Feature weighted ordinal classification) ####
#### Author: Ziyang Ma ####
#### updated on 12/16/2020 ####

## -------- Main function -----------##
FWOC=function(X,Y,W,r,lambda,k,maxiter=100,tol=1e-9){
  ## FWOC (Feature weighted ordinal classification)
  ## Inputs:
  ## X: n x p input matrix
  ## Y: n x K dummy indicator matrix
  ## W: p x p, diagonal matrix denoting weights
  ## r: tuning parameter for weights
  ## lambda: tuning parameter for sparsity
  ## k: number of discriminant vectors
  ## maxiter: number of iterations, default 100
  ## tol: tolerance for convergence
  ## Outputs:
  ## Z: p x k, k discriminant vectors (p x 1)
  
  mats=BW.mat1(X,Y)
  Sw=mats$Sw
  Sb=mats$Sb
  
  p=dim(X)[2]
  A=Sb
  B=r*Sw+(1-r)*W
  Beta=matrix(0,nrow=p,ncol=k)
  V=eigen(A)$vectors[,1:k]
  Z=V ## intialize
  iter=0
  conv=FALSE
  while (iter<maxiter && conv==FALSE){
    Zold=Z
    for (i in 1:p){
      ag= V[i,]-B[i,]%*%Z+B[i,i]%*%Z[i,]
      Z[i,]=ag * max(1-lambda/norm(ag,type="2"),0)/B[i,i]
    }
    if (norm(Zold-Z,type="2")<tol){
      conv=TRUE
    }
    iter=iter+1
  }
  
  return(Z)
}

