# MCMC TVP output
#-----------------

MCMC.TVP.out<-function(fun,burn){
  
  N = fun$N
  TIME = fun$TIME
  P = fun$P
  PA = fun$PA
  
  MCMC.mean.beta<-matrix(NA,TIME,P)
  MCMC.ind<-matrix(NA,TIME,P)
  MCMC.var.beta<-matrix(NA,TIME,P)
  MCMC.mean.alpha<-matrix(NA,TIME,PA)
  MCMC.aux<-matrix(NA,TIME,PA)
  MCMC.var.alpha<-matrix(NA,TIME,PA)
  MCMC.v<-matrix(NA,dim(fun$v)[1],TIME)
  MCMC.var.v<-matrix(NA,dim(fun$v)[1],TIME)
  
  for(i in 1:dim(fun$v)[1]){
    for(t in 1:TIME){
      MCMC.v[i,t] = mean(fun$v[i,t,burn:N])
      MCMC.var.v[i,t] = var(fun$v[i,t,burn:N])
    }
  }
  
  for(p in 1:P){
    for(t in 1:TIME){
      MCMC.mean.beta[t,p] = mean(fun$beta[p,t,burn:N])
      MCMC.var.beta[t,p] = var(fun$beta[p,t,burn:N])
      MCMC.ind[t,p] = mean(fun$gamma_beta[p,t,burn:N])
    }  
  }  
  
  for(p in 1:PA){
    for(t in 1:TIME){
      MCMC.mean.alpha[t,p] = mean(fun$alpha[p,t,burn:N])
      MCMC.var.alpha[t,p] = var(fun$alpha[p,t,burn:N])
      MCMC.aux[t,p] = mean(fun$gamma_alpha[p,t,burn:N])
    }  
  } 
  
  return(list(beta = MCMC.mean.beta,
              alpha = MCMC.mean.alpha,
              ind_beta = MCMC.ind,
              ind_alpha = MCMC.aux,
              var.beta = MCMC.var.beta,
              var.alpha = MCMC.var.alpha,
              v = MCMC.v,
              var.v = MCMC.var.v))
}