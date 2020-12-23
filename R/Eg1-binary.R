devtools::install_github("feizhe/SSHDI")
require(SSHDI)
require(MASS)
library(glmnet);
library(doParallel);
library(foreach);

######### simulate data with binary outcome

n=200; p=300
# n=300; p=500
s0 = c(10,20,30)   #### active set
b0 = rep(0,p)      #### true beta vector
b0[s0] = c(2,2,-2) #### non-zero coefficients

dat = simdat(n,p,xcov = "id",b0=b0,int = -1,family = "binary")
# dat = simdat(n,p,xcov = "ar1",b0=b0,family = "gaussian")
table(dat$y)
dim(dat$x)

m1 = glm(dat$y~ dat$x[,s0], family = "binomial")  ### oracle model
# m1 = glm(dat$y~ dat$x[,s0])
summary(m1)
######### run SSHDI on simulated data

fit1 = SSHDI(dat$x, dat$y, family = "binomial")

###### check outputs
fit1$int                      ####### intercept estimate
##### coefficients and standard error estimates on the active set
cbind(b0, fit1$ss.beta, fit1$sd)[s0,]
hist(-log(fit1$p))            ###### all of derived p-values
fit1$p[s0]                    ###### p-values of the active set signals
fit1$sel.freq[s0]             ###### selection frequency from the SSHDI procedure

