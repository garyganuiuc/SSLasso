EM_lasso <-
function(S,n,p_n,v0,v1,maxiter,p,tau){
  KL1<-NULL
  Q3=NULL
  
  Q1=-Inf
  Q2=100000
  ###initial###
  j=0
  Theta1=100
  W=S
  if(p_n>=n){
    diag(W)=diag(S)+2/n*tau
    Theta=solve(W)
  }else{
    Theta=solve(W)
  }
  
  id=which(upper.tri(Theta)==1)
  
  Theta1=0
  
  #Entropy1=-1000
  P=matrix(0.5,nrow=p_n,ncol=p_n)
  while(j<maxiter&sum(abs(Theta-Theta1))>0.001){
    Q1=-Inf
    tau1=tau
    #update posterior for r
    
    P=1/(1+(v1/v0*exp(-abs(Theta)/(v0)+abs(Theta)/(v1))*(1-p)/p))
    

    Theta1=Theta
    
    
    tau1=tau
    ###update tau#####
    
    for(i in 1:p_n){ 
      Theta1=Theta
      W3=W
      #ensure positive definite
      
      if(tau>0){
        W[i,i]=S[i,i]+2/n*tau
      }
      
      W[-i,i]=W[i,-i]=S[i,-i]+1/(n*v1)*P[i,-i]*sign(Theta[i,-i])+
        1/(n*v0)*(1-P[i,-i])*sign(Theta[i,-i])
      ###update Theta_12#####
      Theta[-i,i]=Theta[i,-i]=
        -Theta[-i,-i]%*%matrix(W[i,-i])/W[i,i]
      
      #  Theta[-i,i]=Theta[i,-i]=v*0.001+Theta[i,-i]
      
      Theta[i,i]=(1-W[i,-i]%*%matrix(Theta[i,-i]))/W[i,i]
      # W[-i,i]=W[i,-i]=
      #  -W[-i,-i]%*%matrix(Theta[i,-i])/Theta[i,i]
      Q2=Qr(Theta,tau,p,P,id,v0,v1,S)
      if(Q2<(Q1)){
        Theta=Theta1
        W=W3
        
        Q2=Q1
      }
      Q1=Q2
      
      
      
    }
    
    
    j=j+1
    # update progress bar
    
  }
  P=1/(1+(v1/v0*exp(-abs(Theta)/(v0)+abs(Theta)/(v1))*(1-p)/p))
  
  return(list(Theta=Theta,P=P,Q=Q2,tau=tau,p=p,W=W))
}
