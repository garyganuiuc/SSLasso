BIC_SSLasso <-
function(Theta,S,P,n){
  p_n=nrow(Theta)
  P1=diag(1,p_n)
  P1[P>0.5]=1
  if(min(eigen(Theta*P1)$values)>0){
    n*(sum(diag(S%*%(Theta*P1)))-log(det(Theta*P1)))+(log(n))*(sum(P>0.5)/2)
  }else{
    Inf
  }
}
