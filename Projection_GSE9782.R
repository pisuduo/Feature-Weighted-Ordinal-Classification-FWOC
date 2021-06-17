#### Application Example for FWOC (GSE9782) ####
#### Author: Ziyang Ma ####

#### Get the GSE9782 Dataset ####
library(Biobase)
library(GEOquery)

#### access data ####
gse9782_data <- getGEO('GSE9782', destdir=".")

expression = gse9782_data$`GSE9782-GPL96_series_matrix.txt.gz`@assayData
expression = expression$exprs  ## matrix of expression data

## treatment label: DEX(control) and PS341
trt_label = gse9782_data[["GSE9782-GPL96_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.1"]]

## response label
response_label=gse9782_data[["GSE9782-GPL96_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.7"]]

## index for use:
all_labels = gse9782_data[["GSE9782-GPL96_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.7"]]
index = intersect(which(trt_label=="treatment = PS341"), which (response_label !="PGx_Response = IE"))

ylab_name= all_labels[index]

ylab_name = droplevels(ylab_name)

## Data matirx n x p
X = expression[,index] ## 22283 x 169
X = t(X)
## create labels
y = rep(0, length(ylab_name))
for ( i in 1:length(ylab_name)){
  if (ylab_name[i]=="PGx_Response = CR"){
    y[i]=1
  }
  if (ylab_name[i]=="PGx_Response = PR"){
    y[i]=2
  }
  if (ylab_name[i]=="PGx_Response = MR"){
    y[i]=3
  }
  if (ylab_name[i]=="PGx_Response = NC"){
    y[i]=4
  }
  if (ylab_name[i]=="PGx_Response = PD"){
    y[i]=5
  }
}

#### ---- Pre-screen the top 500 probes from the orginal dataset ---- ####

## standardization
X_std = scale(X,center=TRUE,scale=TRUE)
pvals = uni.ordilg(X_std,y) ## get the p values for each probe
o = order(pvals,decreasing = FALSE) ## sorted order
ntop=500
feature.ind=o[1:ntop] 
## selected top 500 features
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
fwoc_beta = fwoc(X_sel,Y,W,r=0.8352632,lambda = 0.03,k=2) ## best parameter based on tau and wt_acc
fwoc_beta =ortho(fwoc_beta) ## 

## projections of the original data
proj = X_sel %*% fwoc_beta

## Implementation of penalizedLDA (plda)
library(penalizedLDA)
cv.out = PenalizedLDA.cv(X_sel,y,type="standard",lambdas=c(1e-4,1e-3,1e-2,.1,1.10))
plda_mod = PenalizedLDA(X_sel,y,xte=NULL,lambda=cv.out$bestlambda,K=cv.out$bestK) ## best lambda=0.01, best k=2
res.plda=ortho(plda_mod$discrim)
proj.plda=X_sel%*%res.plda
## get the projection plot Fig 4:

proj_fwoc_df =data.frame(dim1 =proj[,1],dim2=proj[,2],y=y)
proj_plda_df =data.frame(dim1 =proj.plda[,1],dim2=proj.plda[,2],y=y)

library(ggplot2)
p1 = ggplot(data=proj_fwoc_df,aes(x = dim1, y = dim2,color=as.factor(y),shape = as.factor(y)))+
  geom_point()+scale_shape_manual(name = "Class",
                                  labels = c("1", "2", "3",'4','5'),
                                  values = c(15,1,17, 3,19))+
  scale_color_manual(name = "Class",
                     labels = c("1", "2", "3",'4','5'),
                     values = c("#ff6666","#999933","#339966","#3399ff","#cc33cc"))+
  theme_bw()+labs(x="First Dimension",y='Second Dimension',color='Class label')+
  theme(plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))+
  ggtitle("FWOC's projection: 1st vs. 2nd dimensions")+
  theme(legend.position="top")


p2 = ggplot(data=proj_plda_df,aes(x = dim1, y = dim2,color=as.factor(y),shape = as.factor(y)))+
  geom_point()+  scale_shape_manual(name = "Class",
                                    labels = c("1", "2", "3",'4','5'),
                                    values = c(15,1, 17, 3,19))+
  scale_color_manual(name = "Class",
                     labels = c("1", "2", "3",'4','5'),
                     values = c("#ff6666","#999933","#339966","#3399ff","#cc33cc"))+
  theme_bw()+labs(x="First Dimension",y='Second Dimension',color='Class label')+
  theme(plot.title = element_text(hjust=0.5),
        legend.title = element_text(size=12),
        legend.text=element_text(size=12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"))+
  ggtitle("PLDA's projection: 1st vs. 2nd dimensions")+
  theme(legend.position="top")


