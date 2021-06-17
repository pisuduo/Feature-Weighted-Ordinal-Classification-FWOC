### scripts for implementing the method of PLDA, BHGLM, PCRM and FWOC
### prerequisite: run Helper Functions.R and load training and test dataset 
### as matrices: X_tr, y_tr, X_te, y_te

## ------------- implementing PLDA ------------------##
library(penalizedLDA)
library(MASS)
library(caret)
## cross validation
cv.out = PenalizedLDA.cv(X_tr,y_tr,type="standard",lambdas=c(1e-4,1e-3,1e-2,.1,1,10))
## fit model
witten_mod = PenalizedLDA(X_tr,y_tr,xte=X_te,lambda=cv.out$bestlambda,K=cv.out$bestK)
## orthogonal discriminating vectors
res.wit=ortho(witten_mod$discrim) 
## get projections
tr_proj.wit=X_tr%*%res.wit
te_proj.wit=X_te%*%res.wit
## test accuracy
ldamod.wit=lda(tr_proj.wit,grouping=y_tr)
pred.wit=predict(ldamod.wit,te_proj.wit)$class
(acc.wit=mean(pred.wit==y_te))
(tau.wit=cor(as.numeric(pred.wit),y_te,method="kendall")) ## kendall's tau
(wtcost.wit=weightcost(y_te,as.numeric(pred.wit)))           ## weighted cost

(nz.wit=length(unique(c(unlist(apply(res.wit,2,function(x) which(x!=0)))))))

nz.feature.wit=colnames(X_tr)[unique(c(unlist(apply(res.wit,2,function(x) which(x!=0)))))] ## number of non-zero features

##---------- implementing PCRM ---------------------##
library(ordinalgmifs)
Xy_tr=data.frame(X_tr,y_tr)  ## create a dataframe for training set
Xy_te=data.frame(X_te,y_te)  ## create a dataframe for test set

colnames(Xy_te)=colnames(Xy_tr)
penal_cols=names(Xy_tr)[-dim(Xy_tr)[2]]
fit = ordinalgmifs(y_tr ~ 1, x = penal_cols,
                   data = Xy_tr, epsilon = 0.01)

phat <- predict(fit,newx=X_te)
## test accuracy
acc.archer=mean(phat$class==y_te)
wtcost.archer=weightcost(y_te,as.numeric(phat$class))
tau.archer=cor(as.numeric(phat$class),y_te,method="kendall")

nz.archer=length(which(coef(fit)!=0))

##--------implementing bhglm -----------------------##
library(BhGLM)
ytr=as.factor(y_tr)
yte=as.factor(y_te)
Xy_tr=data.frame(X_tr,ytr)  ## create a dataframe for training set
Xy_te=data.frame(X_te,yte)  ## create a dataframe for test set
bpmod=bpolr(ytr~., data=Xy_tr,prior = Student(0,scale=0.1)) #0.1 for result showing
pred=predict(bpmod,X_te)
coef_bh = bpmod$coefficients
### test accuracy

acc.bh=mean(pred==y_te)
tau.bh=cor(as.numeric(pred),y_te,method="kendall") ## kendall's tau
wtcost.bh=weightcost(y_te,as.numeric(pred),1)           ## weighted cost
nz.bh=length(which(coef_bh!=0))

## ---implementing FWOC ----------------------##
J=dim(Y_tr)[2]
max_l=l.max(X_tr,Y_tr,J-1)
lambda_vec=seq(0.01,0.1,by=0.01)
r_vec=seq(0.01,0.99,length=20)
kk=2 ## number of dimension used

### five fold cross validation
cv.res_new=cv_swrlda_gl(X_tr,y_tr,nk=5,r_vec,d=kk,lambda_vec,weight_meth=weightmeth,sel="2",classifier = "KNN")
bestl.swrlda=cv.res_new$best.tau.lambda
bestr.swrlda=cv.res_new$best.tau.r
cvtau=max(cv.res_new$cv.tau)

weight=abs(apply(X_tr,2,cor,y=y_tr,method=weightmeth)) 
W=diag(1-weight,ncol=dim(X_tr)[2],nrow=dim(X_tr)[2])
res.swrlda=S.weight.rLDA.gl(X_tr,Y_tr,W,r=bestr.swrlda,lambda = bestl.swrlda,k=kk)
## discriminating vectors
res.swrlda=ortho(res.swrlda)
## projections
tr_proj=X_tr%*%res.swrlda
te_proj=X_te%*%res.swrlda

### test accuracy
ldamod.swrlda=lda(tr_proj,grouping=y_tr)
pred.swrlda=predict(ldamod.swrlda,te_proj)$class
(acc_lda=mean(pred.swrlda==y_te))
(tau_lda=cor(as.numeric(pred.swrlda),y_te,method="kendall")) ## kendall's tau
wtcost_lda=weightcost(y_te,as.numeric(pred.swrlda),1)           ## weighted cost

(nz.swrlda=length(unique(c(unlist(apply(res.swrlda,2,function(x) which(x!=0)))))))

nz.feature=colnames(X_tr)[which(res.swrlda[,1]!=0)] 

