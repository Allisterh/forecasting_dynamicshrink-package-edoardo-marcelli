# This function performs MAP smoothing in 
# BSTS regression models with stochastic volatility 
# and quarterly seasonality and stochastic trend using 
# Dynamic EMVS with a Particle Smoothing routine
# -----------------------------------------------

DEMVS.Quarterly = function(y,X,N,n0,d0,phi0,phi1,THETA,lambda0,lambda1,delta,
                           sig2u=0.1,sig2d=0.1,sig2tau=0.1,
                           new.X){
  
  require(matlib)
  cat("Starting EMVS for dynamic Shrinkage in Dynamic Linear Model")
  
  m0d=0
  C0d=100
  m0u=0
  C0u=100
  m0tau=0
  C0tau=100
  
  P<-dim(X)[2]
  TIME<-length(y)
  
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P,TIME+1))
  store_v <-matrix(NA,N,TIME+1)
  store_u <-matrix(NA,N,TIME+1)
  store_d <-matrix(NA,N,TIME+1)
  store_tau <-matrix(NA,N,TIME+1)
  
  ols = lm(y ~ X)
  beta0 = ols$coefficients[-1]
  theta<-matrix(beta0,nrow=P,ncol=TIME+1,byrow=F)
  u <- rep(0,TIME+1)
  d <- rep(0,TIME+1)
  tau <- rep(0,TIME+1)
  
  nn<-c()
  dd<-c()
  nn[1]<-n0
  dd[1]<-d0
  v<-c()
  precision<-c()
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  if(missing(new.X)){}else{
    GG.1 = diag(1,P+1+1)
    FF.1 = matrix(c(new.X,1,1),nrow=1)
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
    
    res=y-bigX%*%thetabeta-u[2:(TIME+1)]-tau[2:(TIME+1)]
    if(i %% 100==0)print(res)
    
    for(t in 2:(TIME+1)){
      nn[t]=delta*nn[t-1]+1
      dd[t]=delta*dd[t-1]+res[t-1]^2
    }
    
    precision[TIME+1]=nn[TIME+1]/dd[TIME+1]
    
    v[TIME+1]=1/precision[TIME+1]
    
    for(s in (TIME+1):2){
      
      precision[s-1]=(1-delta)*nn[s-1]/dd[s-1]+delta*precision[s]
      
      v[s-1]=1/precision[s-1]
    }
    
    precision = precision[2:(TIME+1)]
    
    v[1:TIME+1]=0.25
    precision[1:TIME]=4
    
    # COEFFICIENT
    
    for(t in 2:(TIME)){
      
      x = matrix(X[t-1,])
      
      d[t] = ((1/sig2u+2/sig2d)^(-1))*((1/sig2u)*(u[t+1]-u[t])+(1/sig2d)*(d[t-1]+d[t+1]))
      
      u[t] = ((sig2u*v[t])/(2*v[t]+sig2u))*(1/sig2u*(u[t+1]+u[t-1]+d[t-1]-d[t])+1/v[t]*(y[t-1]-theta[,t]%*%x-tau[t]))
      
      D = diag(gamma[t,]/lambda1+(1-gamma[t,])/lambda0+
                 (phi1^2)*gamma[t+1,]/lambda1)
      
      X1 = x%*%t(x)
      
      Sigma = precision[t-1]*X1+D
      
      invD = diag(1/diag(D))
      
      Den = 1 + precision[t-1]*t(x)%*%invD%*%x
      
      invDen = 1/Den[1,1]
      
      invSigma = invD - precision[t-1]*invD%*%(X1*invDen)%*%invD
      
      nu=precision[t-1]*(y[t-1]-u[t]-tau[t])*x+phi1/lambda1*theta[,t-1]*gamma[t,]+(phi1/lambda1*theta[,t+1]*gamma[t+1,])
      
      theta[,t] = invSigma%*%nu
      
    }
    
    x = matrix(X[TIME,])
    
    d[TIME+1] = ((1/sig2d)^(-1))*((1/sig2d)*(d[TIME]))
    
    u[TIME+1] = ((sig2u*v[TIME])/(1*v[TIME]+sig2u))*(1/sig2u*(u[TIME]+d[TIME])+1/v[TIME]*(y[TIME]-theta[,TIME+1]%*%x-tau[TIME+1]))
    
    D = diag(gamma[TIME+1,]/lambda1+(1-gamma[TIME+1,])/lambda0)
    
    X1 = x%*%t(x)
    
    Sigma=precision[TIME]*x%*%t(x)+D
    
    nu=precision[TIME]*(y[TIME]-u[TIME+1]-tau[TIME+1])*x+phi1/lambda1*theta[,TIME]*gamma[TIME+1,]
    
    invD = diag(1/diag(D))
    
    Den = 1 + precision[t]*t(x)%*%invD%*%x
    
    invDen = 1/Den[1,1]
    
    invSigma = invD - precision[TIME]*invD%*%(X1*invDen)%*%invD
    
    theta[,TIME+1]= invSigma%*%nu
    
    Sigma = diag((1-phi1^2)*gamma[1,]/lambda1+(1-gamma[1,])/lambda0+(phi1^2)*gamma[2,]/lambda1)
    
    invSigma = diag(1/diag(Sigma))
    
    theta[,1] = phi1/lambda1*invSigma%*%theta[,2]*gamma[2,]
    
    d[1]=((1/C0d+1/sig2u+1/sig2d)^(-1))*(m0d/C0d+1/sig2u*(u[2]-u[1])+1/sig2d*d[2])
    
    u[1]=((1/C0u+1/sig2u)^(-1))*(m0u/C0u+u[2]/sig2u-d[1]/sig2u)
    
    # Seasonality
    
    tau[2] = ((4/sig2tau+1/v[2])^(-1))*(1/sig2tau*(-3*tau[1]-3*tau[3]-2*tau[4]-tau[5])
                                        +1/v[2]*(y[1]-theta[,2]%*%matrix(X[1,])-u[2]))
    
    tau[3] = ((4/sig2tau+1/v[3])^(-1))*(1/sig2tau*(-3*tau[2]-2*tau[1]-3*tau[4]-2*tau[5]-tau[6])
                                        +1/v[3]*(y[1]-theta[,3]%*%matrix(X[2,])-u[3]))
    
    tau[4] = ((4/sig2tau+1/v[4])^(-1))*(1/sig2tau*(-3*tau[3]-2*tau[2]-tau[1]-3*tau[5]-2*tau[6]-tau[7])
                                        +1/v[4]*(y[1]-theta[,4]%*%matrix(X[4,])-u[4]))
    
    
    tau[TIME+1] = ((1/sig2tau+1/v[TIME+1])^(-1))*(1/sig2tau*(-tau[TIME]-tau[TIME-1]-tau[TIME-2])
                                                  +1/v[TIME+1]*(y[TIME]-theta[,TIME+1]%*%matrix(X[TIME,])-u[TIME+1]))
    
    tau[TIME] = ((2/sig2tau+1/v[TIME])^(-1))*(1/sig2tau*(-2*tau[TIME-1]-2*tau[TIME-2]-tau[TIME-3]-tau[TIME+1])
                                              +1/v[TIME]*(y[TIME-1]-theta[,TIME]%*%matrix(X[TIME-1,])-u[TIME]))
    
    tau[TIME-1] = ((4/sig2tau+1/v[TIME-1])^(-1))*(1/sig2tau*(-3*tau[TIME-2]-2*tau[TIME-3]-tau[TIME-4]-2*tau[TIME]-tau[TIME+1])
                                                  +1/v[TIME-1]*(y[TIME-2]-theta[,TIME-1]%*%matrix(X[TIME-2,])-u[TIME-1]))
    
    tau[TIME-2] = ((4/sig2tau+1/v[TIME-2])^(-1))*(1/sig2tau*(-3*tau[TIME-3]-2*tau[TIME-4]-tau[TIME-5]-3*tau[TIME-1]-2*tau[TIME]-tau[TIME+1])
                                                  +1/v[TIME-2]*(y[TIME-3]-theta[,TIME-2]%*%matrix(X[TIME-3,])-u[TIME-2]))
    
    for(t in 5:(TIME-3)){
      
      tau[t] = ((4/sig2tau+1/v[t])^(-1))*(1/sig2tau*(-3*tau[t-1]-2*tau[t-2]-tau[t-3]-3*tau[t+1]-2*tau[t+2]-tau[t+3])
                                          +1/v[t]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    }
    
    tau[1] = ((3/sig2tau+1/C0tau)^(-1))*(m0tau/C0tau+1/sig2tau*(-3*tau[2]-2*tau[3]-tau[4]))
    
    # STORE SAMPLES
    
    store_gamma[i,,] = t(gamma) 
    store_theta[i,,] = theta 
    store_v[i,] = v
    store_d[i,] = d
    store_u[i,] = u
    store_tau[i,] = tau
  }
  
  if(missing(new.X)){ }else{
    
    state_vec = matrix(c(theta[,TIME+1],u[TIME+1],tau[TIME+1]),ncol=1)
    a = GG.1%*%state_vec
    f = FF.1%*%a
    
  }
  
  cat("MCMC takes")
  print(proc.time() - ptm)
  
  cat("")
  
  return(list(beta=store_theta,gamma=store_gamma,v=store_v,u=store_u,d=store_d,tau = store_tau,
              fc.f = f))
  
}
