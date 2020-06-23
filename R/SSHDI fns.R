#
# Author: Zhe Fei
# Date: Nov. 12, 2019
#
# Estimation and inference for high-dimensional linear model with SPARES
# "Inference for High-dimensional Linear Model: A Selection-assisted Partial Regression and Smoothing Approach"
#
#
library(MASS);
library(glmnet);
library(doParallel);
library(foreach);
#######

extr<-function(fit1=list(),k,pos){
  temp1<-fit1[[k]][[pos]]
}

est.var<- function(beta,Ycount, unbiased = TRUE){
  n = dim(Ycount)[1]
  B = dim(Ycount)[2]
  n1 = sum(Ycount[,1])

  var1 <- (n/(n-n1))^2 * sum(cov(beta,t(Ycount))^2)
  var2 <- ifelse(unbiased, var1 - n1/(n-n1)*n/B*var(beta), var1)
  return(var2)
}

glmj<-function(j,yvec,xmat,samp1,set,family = "gaussian"){
  fitj<-glm(yvec[samp1]~xmat[samp1,union(j,set)],family = family)
  as.numeric(coef(fitj)[2])
}



SPARE<-function(xmat,yvec,family = "gaussian", prop = 1/2,p0=NA ,lam0 = NA){

  n = dim(xmat)[1]
  p = dim(xmat)[2]
  if ( n != length(yvec) ) print("Lengths differ.")
  n1 = n*prop
  p0 = ifelse(is.na(p0), sqrt(p), p0)

  samp1<-sort(sample(1:n,n/2))
  samp2<-(1:n)[-samp1]

  l2<-glmnet(xmat[samp2,],yvec[samp2], family = family )
  idx <- ifelse(is.na(lam0), max(which(l2$df<p0)), which(l2$lambda<lam0)[1])

  coef1<-l2$beta[,idx]
  set<-(1:p)[coef1!=0]

  beta_hat<-rep(0,p)

  coef2<-coef(glm(yvec[samp1]~xmat[samp1,set],family = family ))
  beta_hat[set]<-coef2[-1]
  beta_hat[-set]<-sapply((1:p)[-set],glmj,yvec=yvec,xmat=xmat,samp1=samp1,set=set,family = family)

  yicount<-rep(0,n)
  yicount[samp1] <- 1

  returnList<-list("int" = coef2[1],"beta.hat"=beta_hat,"sel.set"=set,"boot.ct"=yicount)
  return(returnList)
}


SSHDI<-function(xmat,yvec,family = "gaussian",B=500){

  glmj<-function(j,yvec,xmat,samp1,set,family){
    fitj<-glm(yvec[samp1]~xmat[samp1,union(j,set)],family = family)
    as.numeric(coef(fitj)[2])
  }

  SPARE<-function(xmat,yvec,family, prop = 1/2,p0=NA ,lam0 = NA){

    n = dim(xmat)[1]
    p = dim(xmat)[2]
    n1 = n*prop
    p0 = ifelse(is.na(p0), sqrt(p), p0)

    samp1<-sort(sample(1:n,n/2))
    samp2<-(1:n)[-samp1]

    l2<-glmnet(xmat[samp2,],yvec[samp2], family = family )
    idx <- ifelse(is.na(lam0), max(which(l2$df<p0)), which(l2$lambda<lam0)[1])

    coef1<-l2$beta[,idx]
    set<-(1:p)[coef1!=0]

    beta_hat<-rep(0,p)

    coef2<-coef(glm(yvec[samp1]~xmat[samp1,set],family = family ))
    beta_hat[set]<-coef2[-1]
    beta_hat[-set]<-sapply((1:p)[-set],glmj,yvec=yvec,xmat=xmat,samp1=samp1,set=set,family = family)

    yicount<-rep(0,n)
    yicount[samp1] <- 1

    returnList<-list("int" = coef2[1],"beta.hat"=beta_hat,"sel.set"=set,"samp.ct"=yicount)
    return(returnList)
  }

  n<-dim(xmat)[1]
  p<-dim(xmat)[2]

  l0<-cv.glmnet(xmat,yvec,family = family)
  lam0<-l0$lambda.1se

  print("running...")
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  fit1<- foreach(i = 1:B,.packages=c("MASS","glmnet")) %dopar% SPARE(xmat,yvec,family = family,p0 = 9)
  stopCluster(cl)
  print("done")
  INT<-sapply(1:B,extr,fit1=fit1,pos=1)
  BETA<-sapply(1:B,extr,fit1=fit1,pos=2)
  SET<-lapply(1:B,extr,fit1=fit1,pos=3)
  Ycount<-sapply(1:B,extr,fit1=fit1,pos=4)

  betam<-apply(BETA,1, median)
  vars<-apply(BETA,1,est.var,Ycount=(Ycount))
  sds <- sqrt(vars)
  pvs<-2*(1-pnorm(abs(betam)/sds))
  nzero <- (1:p)[pvs<0.05/p & !is.na(pvs)]
  int<-median(INT)

  temptab<-table(unlist(SET))
  sel_freq<-rep(0,p)
  sel_freq[as.numeric(names(temptab))]<-temptab/B

  returnList<-list("int"=int,"ss.beta"=betam,"sd"=sds,"p"=pvs,"sel.freq"=sel_freq,
                   "nzero"=nzero)
  return(returnList)
}
