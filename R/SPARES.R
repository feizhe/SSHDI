#
# Author: Zhe Fei
# Date: March 11, 2017
#
# Estimation and inference for high-dimensional linear model with SPARES
# "Inference for High-dimensional Linear Model: A Selection-assisted Partial Regression and Smoothing Approach"
#
#
# library(Matrix);
# library(expm);
# library(flare);
# library(Hmisc)
# library(lars)
# library(SIS)
#
library(MASS);
library(glmnet);
library(doParallel);
library(foreach);
#######

extr<-function(fit1,k,pos){
  temp1<-fit1[[k]][[pos]]
}

sd_smooth<-function(beta,Ycount){
  return(sqrt(sum(cov(beta,t(Ycount))^2)))
}

sdu<-function(beta,Ycount,nloop,n){
  sd2<-sum(cov(beta,t(Ycount))^2)
  return(sqrt(sd2-n/2/nloop*var(beta)))
}

SPARES<-function(xmat,yvec,selection="lasso",nloop=500){
  PaRes<-function(x1,hs,m){
    t(x1)%*%(diag(rep(1,m))-hs)%*%x1
  }
  
  SPARE<-function(xmat,yvec,selection="lasso",n,p,lam0){
    
    samp1<-sample(1:n,n/2,replace=T)
    tab1<-table(samp1)
    samp2<-(1:n)[-as.numeric(names(tab1))]
    n2<-length(samp2)
    
    l2<-glmnet(xmat[samp2,],yvec[samp2])
    lam2<-(1:length(l2$lambda))[l2$lambda<lam0][1]
    coef2<-l2$beta[,lam2]
    set<-(1:p)[coef2!=0]
    
    beta_hat<-rep(0,p)
    xs<-xmat[samp1,set]
    pxs<-solve(t(xs)%*%xs)
    hs<-xs%*%pxs%*%t(xs)
    beta_hat[set]<-pxs%*%t(xs)%*%yvec[samp1]
    
    x_s<-xmat[samp1,-set]
    h_s<-apply(x_s,2,PaRes,hs=hs,m=n/2)
    beta_hat[-set]<-(t(x_s)/h_s)%*%(diag(rep(1,n/2))-hs)%*%yvec[samp1]
    yicount<-rep(0,n)
    yicount[-samp2]<-tab1
    
    returnList<-list("beta.hat"=beta_hat,"sel.set"=set,"boot.ct"=yicount)
    return(returnList)
  }
  
  n<-dim(xmat)[1]
  p<-dim(xmat)[2]
  
  l0<-cv.glmnet(xmat,yvec)
  lam0<-l0$lambda.1se
  
  y0<-mean(yvec)
  yvec1<-yvec-y0
  xmat1<-scale(xmat)
  sdx<-apply(xmat,2,sd)
  
  print("running...")
  registerDoParallel(cores=4)
  fit1<- foreach(i = 1:nloop,.packages=c("MASS","glmnet")) %dopar% SPARE(xmat1,yvec1,selection="lasso",n,p,lam0)
  print("done")
  BETA<-sapply(1:nloop,extr,fit1=fit1,pos=1)
  SET<-lapply(1:nloop,extr,fit1=fit1,pos=2)
  Ycount<-sapply(1:nloop,extr,fit1=fit1,pos=3)
  
  betam<-rowMeans(BETA)
  sds<-apply(BETA,1,sdu,Ycount=(Ycount),nloop=nloop,n=n)
  pvs<-2*(1-pnorm(abs(betam)/sds))
  nzero<-(1:p)[pvs<0.05/p]
  int<-y0-sum(betam[nzero]*colMeans(xmat[,nzero])/sdx[nzero])
  
  temptab<-table(unlist(SET))
  sel_freq<-rep(0,p)
  sel_freq[as.numeric(names(temptab))]<-temptab/nloop
  
  returnList<-list("sm.beta"=betam/sdx,"sd"=sds/sdx,"p"=pvs,"sel.freq"=sel_freq,
                   "int"=int,"nzero"=nzero)
  return(returnList)
}