#---------------------------------------------------------------------------#
#                                                                           #
# Sparse Time-Series Generator                                              #
#                                                                           #
#---------------------------------------------------------------------------#


sparse.data.sim<-function(TIME,P,FP,phi0,phi1,lambda1,v,seed){
  
  if(missing(seed)){
    cat("No seed has been setted")
  }else{
    set.seed(seed)
  }
  b<-matrix(NA,nrow=TIME,ncol=P)
  zerob<-matrix(NA,nrow=TIME,ncol=FP)
  c<-matrix(NA,nrow=TIME,ncol=P)
  X<-matrix(NA,nrow=TIME,ncol=P)
  true_ind<-matrix(NA,nrow=TIME,ncol=(P+FP))
  y<-c()
  
  for(p in 1:P){
    i = rbinom(1,1,0.5)
    if(i==1){
      b[1,p] = runif(1,2,3)
    }else{
      b[1,p] = runif(1,-3,-2)
    }
  }
  
  for(p in 1:P){
    X[1,p]<-rnorm(1,0,1)
  }
  
  y[1]<-t(X[1,])%*%b[1,]*2+rnorm(1,0,sqrt(v))
  
  for(t in 2:TIME){
    
    for(p in 1:P){
      c[t,p]<-phi0+phi1*(b[t-1,p]-phi0)+rnorm(1,0,sqrt(lambda1))
      if(c[t,p]>-0.5 & c[t,p]<0.5){b[t,p]=0}else{b[t,p]=c[t,p]}
    }
    
    for(p in 1:P){
      X[t,p]<-rnorm(1,0,1)
    }
    
    y[t]<-t(X[t,])%*%b[t,]*2+rnorm(1,0,sqrt(v))
    
  }
  
  if(missing(FP)){}else{
    X.false<-matrix(NA,nrow=TIME,ncol=FP)
    
    for(t in 1:TIME){  
      for(p in 1:FP){
        X.false[t,p]<-rnorm(1,0,1)
        zerob[t,p] = 0
      }
    }
    
    X<-cbind(X,X.false)
  }
  
  true_b = cbind(b*2,zerob)
  for(p in 1:(P+FP)){
    for(t in 1:TIME){
      if(true_b[t,p]==0){true_ind[t,p]=0}else{true_ind[t,p]=1}
    }
  }
  
  
  return(list(y=y,X=X,b=(b*2),true_ind=true_ind,true_b=true_b))
}



