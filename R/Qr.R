Qr <-
function(Theta,tau,p,P,id,v0,v1,S){
  n=nrow(Theta)
  n/2*log(det(Theta))-n/2*sum(diag(S%*%Theta))+
    (-tau*sum(diag(Theta))+sum((-log(2*v1)-abs(Theta[id])/v1)*P[id]+
                                (-log(2*v0)-abs(Theta[id])/v0)*(1-P[id]))+
    sum( P[id]*log(p)+(1-P[id])*log(1-p)))
}
