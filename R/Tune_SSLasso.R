Tune_SSLasso <-
function(v0,tau,S,n,p_n,p){
pb <- txtProgressBar(min = 0, max = length(v0)*10*length(tau), style = 3)
bic=NULL
for(m in 1:length(p)){
  for(i in 1:length(v0)){
    v1=seq(v0[i]+0.1,5*v0[i],0.5*v0[i])
      for(j in 1:length(v1)){
        for(k in 1:length(tau)){
          ##Bayes EM####
          #v1=sqrt(n/(log(p_n)))/5
          #v0=sqrst(n/(p_n*log(p_n)))/5
          
          #tau=sqrt((log(p_n)/n))
          w=1
          l=1
          maxiter=30
          result1<-EM_lasso(S,n,p_n,v0[i],v1[j],maxiter,p[m],tau[k])
          
          bic=rbind(bic,list(v0=v0[i],v1=v1[j],tau=tau[k],p=p[m],BIC=BIC_SSLasso(result1$Theta,S,result1$P,n)))
          
          setTxtProgressBar(pb, (m-1)*length(v0)*length(v1)*length(tau)+(i-1)*length(v1)*length(tau)+(j-1)*length(tau)+k)
        
      }
    } 
  }
  
}

close(pb)
return(bic[which.min(bic[,5]),])
}
