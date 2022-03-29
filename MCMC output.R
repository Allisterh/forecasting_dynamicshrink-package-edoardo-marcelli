# MCMC output
#-------------

MCMC.out<-function(fun,burn){
  
  N<-dim(fun$theta)[1]
  PS<-dim(fun$theta)[2]
  P<-dim(fun$gamma)[2]
  TIME<-dim(fun$fit.reg)[2]
  
  MCMC.mean<-matrix(NA,TIME,PS)
  #MCMC.C<-matrix(NA,TIME,PS)
  MCMC.indicator<-matrix(NA,TIME,P)
  MCMC.var<-matrix(NA,TIME,PS)
  MCMC.f<-matrix(NA,TIME)
  MCMC.v<-matrix(NA,TIME)
  MCMC.var.v<-matrix(NA,TIME)
  #MCMC.q<-matrix(NA,TIME)
  MCMC.varf<-matrix(NA,TIME)
  
  
  for(t in 1:TIME){
    MCMC.f[t]<-mean(fun$fit.reg[burn:N,t])
    #MCMC.q[t]<-mean(fun$q[burn:N,t])
    MCMC.v[t]<-mean(fun$v[burn:N,t])
    MCMC.var.v[t]<-var(fun$v[burn:N,t])
    MCMC.varf[t]<-var(fun$fit.reg[burn:N,t])
    
    for(p in 1:P){
      MCMC.indicator[t,p]<-mean(fun$gamma[burn:N,p,t+1])
    }
    
    for(p in 1:PS){
      MCMC.mean[t,p]<-mean(fun$theta[burn:N,p,t+1])
      MCMC.var[t,p]<-var(fun$theta[burn:N,p,t+1])
    } 
    
    
    
  }
  
  return(list(mean=MCMC.mean,fit.reg=MCMC.f,ind=MCMC.indicator,
              var=MCMC.var,v=MCMC.v,var.v = MCMC.var.v, var.fit.reg=MCMC.varf))
}
