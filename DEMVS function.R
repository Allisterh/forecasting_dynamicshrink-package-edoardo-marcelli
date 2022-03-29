# Dynamic EMVS for TVp regression model with discount factor model for variances
#--------------------------------------------------------------------------------

DEMVS<-function(y,X,N,n0,d0,phi0,phi1,THETA,lambda0,lambda1,delta,new.X){
  
  require(matlib)
  cat("Starting EMVS for dynamic Shrinkage in Dynamic Linear Model")
  
  # Preliminaries
  #-----------------------------------------
  
  P<-dim(X)[2]
  TIME<-length(y)
  
  ols = lm(y ~ X)
  beta0 = ols$coefficients[-1]
  beta0[beta0>10 | beta0<(-10)]=10
  beta0[is.na(beta0)]=0
  theta<-matrix(beta0,nrow=P,ncol=TIME+1)
  n<-c()
  d<-c()
  n[1]<-n0
  d[1]<-d0
  v<-c()
  precision<-c()
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  # Storing Objects
  #------------------------------
  
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P,TIME+1))
  store_v <-matrix(NA,N,TIME)
  
  # EM
  #---------------------------------
  
  if(missing(new.X)){}else{
    FF.1 = matrix(c(new.X),nrow=1)
  } 
  
  ptm <- proc.time()
  
  for(i in 1:N){
    
    if(i %% 100==0) cat("\n", i, "loops...") 
    if(i == N) cat("Done!\n") 
    
    # INDICATOR
    
    beta<-matrix(theta,ncol=1,byrow=T)
    thetabeta<-matrix(theta[,2:(TIME+1)],ncol=1,byrow=T)
    
    thetaind = THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                    +(1-THETA)*dnorm(beta[1:((TIME+1)*P-P),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(beta[1:((TIME+1)*P-P),1]-phi0)
    
    pstar = thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))/(thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))
                                                                         +(1-thetaind)*dnorm(beta[(P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
    
    pstar0 = thetaind[1:P]
    
    gamma = matrix(c(pstar0,pstar),ncol=P,byrow=T)
    
    # VOLATILITY
    
    res=y-bigX%*%thetabeta
    
    for(t in 2:(TIME+1)){
      n[t]=delta*n[t-1]+1
      d[t]=delta*d[t-1]+res[t-1]^2
    }
    
    precision[TIME+1]=n[TIME+1]/d[TIME+1]
    
    v[TIME+1]=1/precision[TIME+1]
    
    for(s in (TIME+1):2){
      
      precision[s-1]=(1-delta)*n[s-1]/d[s-1]+delta*precision[s]
      
      v[s-1]=1/precision[s-1]
    }
    
    precision = precision[2:(TIME+1)]
    if(i<10){
      precision = rep(4,TIME)
    }else{}
    
    # COEFFICIENT
    
    for(t in 1:(TIME-1)){
      
      x = matrix(X[t,])
      
      D = diag(gamma[t+1,]/lambda1+(1-gamma[t+1,])/lambda0+
                 (phi1^2)*gamma[t+2,]/lambda1)
      
      X1 = x%*%t(x)
      
      Sigma = precision[t]*X1+D
      
      invD = diag(1/diag(D))
      
      Den = 1 + precision[t]*t(x)%*%invD%*%x
      
      invDen = 1/Den[1,1]
      
      invSigma = invD - precision[t]*invD%*%(X1*invDen)%*%invD
      
      u=precision[t]*y[t]*x+phi1/lambda1*theta[,t]*gamma[t+1,]+(phi1/lambda1*theta[,t+2]*gamma[t+2,])
      
      theta[,t+1] = invSigma%*%u
      
    }
    
    D = diag(gamma[TIME+1,]/lambda1+(1-gamma[TIME+1,])/lambda0)
    
    x = matrix(X[TIME,])
    
    X1 = x%*%t(x)
    
    Sigma=precision[TIME]*x%*%t(x)+D
    
    u=precision[TIME]*y[TIME]*x+phi1/lambda1*theta[,TIME]*gamma[TIME+1,]
    
    invD = diag(1/diag(D))
    
    Den = 1 + precision[TIME]*t(x)%*%invD%*%x
    
    invDen = 1/Den[1,1]
    
    invSigma = invD - precision[TIME]*invD%*%(X1*invDen)%*%invD
    
    theta[,TIME+1]= invSigma%*%u
    
    Sigma = diag((1-phi1^2)*gamma[1,]/lambda1+(1-gamma[1,])/lambda0+(phi1^2)*gamma[2,]/lambda1)
    
    invSigma = diag(1/diag(Sigma))
    
    theta[,1] = phi1/lambda1*invSigma%*%theta[,2]*gamma[2,]
    
    # STORE SAMPLES
    
    store_gamma[i,,] = t(gamma) 
    store_theta[i,,] = theta 
    store_v[i,] = v[2:length(v)]
    
  }
  
  if(missing(new.X)){f = NULL }else{
    
    GG.1 = diag(phi1*gamma[dim(gamma)[1],],P)
    state_vec = matrix(c(theta[,TIME+1]),ncol=1)
    a = GG.1%*%state_vec
    f = FF.1%*%a
    
  }
  
  cat("MCMC takes")
  print(proc.time() - ptm)
  
  cat("For further details please consult Rockova & McAlinn 2021")
  
  return(list(beta=store_theta,fc.f = f,gamma=store_gamma,v=store_v))
  
}