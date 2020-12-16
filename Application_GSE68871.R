#### Application Example for FWOC (GSE68871) ####
#### Author: Ziyang Ma ####
#### Updated on 12/16/2020 ####

#### Get the GSE68871 Dataset ####
library(Biobase)
library(GEOquery)

#### access data ####
gse68871_data <- getGEO('GSE68871', destdir=".")
exprs_X=gse68871_data[["GSE68871_series_matrix.txt.gz"]]@assayData[["exprs"]]
exprs_label=gse68871_data[["GSE68871_series_matrix.txt.gz"]]@phenoData@data[["response to vtd therapy:ch1"]]
## create labels
y=rep(0,length(exprs_label))
for(i in 1:length(y)){
  if (exprs_label[i]=="CR"){
    y[i]=1
  }
  if (exprs_label[i]=="nCR"){
    y[i]=2
  }
  if (exprs_label[i]=="VGPR"){
    y[i]=3
  }
  if (exprs_label[i]=="PR"){
    y[i]=4
  }
  if (exprs_label[i]=="SD"){
    y[i]=5
  }
}


## Data matirx n x p
X=t(exprs_X) # 118 x 54675


#### ---- Pre-screen the top 1000 probes from the orginal dataset ---- ####

## standardization
X_std = scale(X,center=TRUE,scale=TRUE)
pvals = uni.ordilg(X_std,y) ## get the p values for each probe
o = order(pvals,decreasing = FALSE) ## sorted order
ntop=1000
feature.ind=o[1:ntop] 
## selected top 1000 features
X_sel = X_std[,feature.ind] 
colnames(X_sel) = sub('.','',colnames(X_sel))

## n x K dummy indicator matrix
Y = create_dummyY(y) 
J=dim(Y)[2]
## get the maximum of lambda 
max_l=l.max(X_sel,Y,J-1) 
lambda_vec=seq(0.01,max_l,by=0.01)
r_vec=seq(0.01,0.99,length=20)
kk=2

### 5-fold cross validation
## cv.res_new=cv_fwoc(X_sel,y,nk=5,r_vec,d=kk,lambda_vec,weight_meth=weightmeth,sel="2",classifier = "KNN")

### calculate feature weight matrix ####
weightmeth="kendall"
weight=abs(apply(X_sel,2,cor,y=y,method=weightmeth)) 
W=diag(1-weight,ncol=dim(X_sel)[2],nrow=dim(X_sel)[2])

### get the discriminating vectors ###
fwoc_beta = fwoc(X_sel,Y,W,r=0.319473684210526,lambda = 0.03,k=2) ## best parameter based on tau and wt_acc
fwoc_beta =ortho(fwoc_beta) ## 

## projections of the original data
proj = X_sel %*% fwoc_beta
