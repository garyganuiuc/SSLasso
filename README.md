# BayesEM
R package (Sparse Precision Matrix Estimation with Bayesian Regularization)

##Maintainer
Lingrui(Gary) Gan, University of Illinois at Urbana-Champaign

##Description
Bayesian Model for learning the estimating sparse Precision Matrix with Statistical Guarantee.

##Example
###Model Set up
```{r}
p_n=p=50 #number of variables
n=100 #number of observations
C=toeplitz(c(1,0.5,rep(0,p_n-2)))
Sigma=solve(C)

Sigma<-solve(C)
Y<-mvrnorm(n,rep(0,p_n),Sigma)
S<-cov(Y)  #sample covariance
```

###Tuning and Estimate
```{r}
v0=c(0.1,0.13)

tau=c(0.04,0.07,0.1,0.4)*n/2

Sigma<-solve(C)
Y<-mvrnorm(n,rep(0,p_n),Sigma)
S<-cov(Y)  #sample covariance
p=0.5
Tune=Tune_SSLasso(v0,tau,S,n,p_n,p)
maxiter=30
v0_t=Tune$v0
v1_t=Tune$v1
tau_t=Tune$tau
result1<-EM_lasso(S,n,p_n,v0_t,v1_t,maxiter,p,tau_t)
```


##Reference
Manuscript  
