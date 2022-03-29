# This function performs MAP smoothing in 
# BSTS regression models with stochastic volatility 
# and monthly seasonality and stochastic trend using 
# Dynamic EMVS with a Particle Smoothing routine
# -----------------------------------------------

DEMVS.Monthly<-function(y,X,N,phi0,phi1,THETA,lambda0,lambda1,
                        sig2u=0.1,sig2d=0.1,sig2tau=0.1,
                        new.X){
  
  require(matlib)
  cat("Starting Dynamic EMVS for Bayesian Structural Time Series models with monthly seasonality")
  
  m0d=0
  C0d=100
  m0u=0
  C0u=100
  m0tau=0
  C0tau=100
  
  # Particle Filtering
  #--------------------
  
  BPF_SV<-function(data,N,alpha0,alpha1,sigzeta,r){
    if(missing(r)){r=2}else{}
    xf = matrix(NA,nrow=length(data),ncol=N)
    ws = matrix(NA,nrow=length(data),ncol=N)
    x = rnorm(N,alpha0,sigzeta/sqrt(1-alpha1^2))
    w = rep(1/N,N)
    for(t in 1:length(data)){
      
      x = alpha0+alpha1*(x-alpha0)+rnorm(N,0,sigzeta)
      w1 = w*dnorm(data[t],0,exp(x/2))
      w = w1/sum(w1)
      ESS = 1/sum(w^2)
      if(ESS<N/r){
        index<-sample(N,size=N,replace=T,prob=w)
        x<-x[index]
        w<-rep(1/N,N)
      }else{}
      xf[t,] = x
      ws[t,] = w
    }
    xs = xf
    for(t in (length(data)-1):1){
      ws[t,] = dnorm(xs[t+1,],alpha0+alpha1*(xf[t,]-alpha0),sigzeta)
      index<-sample(N,size=N,replace=T,prob=ws[t,])
      xs[t,]<-xf[t,index]
    }
    xhat = sapply(1:TIME, function(i) weighted.mean(xs[i,],ws[i,]))
    return(list(xhat=xhat,xs=xs,ws=ws))
  }
  
  # Preliminaries
  #---------------
  
  P<-dim(X)[2]
  TIME<-length(y)
  
  ols = lm(y ~ X)
  beta0 = ols$coefficients[-1]
  theta<-matrix(beta0,nrow=P,ncol=TIME+1,byrow=F)
  u <- rep(0,TIME+1)
  d <- rep(0,TIME+1)
  tau <- rep(0,TIME+1)
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  if(missing(new.X)){}else{
    FF.1 = matrix(c(new.X,1,1),nrow=1)
  } 
  
  # Storing objects
  #------------------
  
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P,TIME+1))
  store_v <-matrix(NA,N,TIME)
  store_u <-matrix(NA,N,TIME+1)
  store_d <-matrix(NA,N,TIME+1)
  store_tau <-matrix(NA,N,TIME+1)
  store_reg <-matrix(NA,N,TIME)
  
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
    reg=bigX%*%thetabeta-u[2:(TIME+1)]
    
    if(i<10){
      v=rep(0.25,TIME)
      precision=rep(4,TIME)
    }else{
      
      out.pf = BPF_SV(res,N=1000,alpha0=-2,alpha1=0.9,sigzeta=0.1)
      h = out.pf$xhat
      v = exp(h)
      precision = 1/v
      
    }
    
    # COEFFICIENT
    
    for(t in 2:(TIME)){
      
      x = matrix(X[t-1,])
      
      d[t] = ((1/sig2u+2/sig2d)^(-1))*((1/sig2u)*(u[t+1]-u[t])+(1/sig2d)*(d[t-1]+d[t+1]))
      
      u[t] = ((sig2u*v[t-1])/(2*v[t-1]+sig2u))*(1/sig2u*(u[t+1]+u[t-1]+d[t-1]-d[t])+1/v[t-1]*(y[t-1]-theta[,t]%*%x-tau[t]))
      
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
    
    Den = 1 + precision[TIME]*t(x)%*%invD%*%x
    
    invDen = 1/Den[1,1]
    
    invSigma = invD - precision[TIME]*invD%*%(X1*invDen)%*%invD
    
    theta[,TIME+1]= invSigma%*%nu
    
    Sigma = diag((1-phi1^2)*gamma[1,]/lambda1+(1-gamma[1,])/lambda0+(phi1^2)*gamma[2,]/lambda1)
    
    invSigma = diag(1/diag(Sigma))
    
    theta[,1] = phi1/lambda1*invSigma%*%theta[,2]*gamma[2,]
    
    d[1]=((1/C0d+1/sig2u+1/sig2d)^(-1))*(m0d/C0d+1/sig2u*(u[2]-u[1])+1/sig2d*d[2])
    
    u[1]=((1/C0u+1/sig2u)^(-1))*(m0u/C0u+u[2]/sig2u-d[1]/sig2u)
    
    # Seasonality
    
    t=TIME+1
    tau[t] = ((1/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t-1]-tau[t-2]-tau[t-3]-tau[t-4]-tau[t-5]
                                                     -tau[t-6]-tau[t-7]-tau[t-8]-tau[t-9]-tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    t=TIME
    tau[t] = ((2/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+1]-2*tau[t-1]-2*tau[t-2]-2*tau[t-3]-2*tau[t-4]-2*tau[t-5]
                                                     -2*tau[t-6]-2*tau[t-7]-2*tau[t-8]-2*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    t = TIME-1
    tau[t] = ((3/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+2]-2*tau[t+1]-3*tau[t-1]-3*tau[t-2]-3*tau[t-3]-3*tau[t-4]-3*tau[t-5]
                                                     -3*tau[t-6]-3*tau[t-7]-3*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    t=TIME-2
    tau[t] = ((4/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+3]-2*tau[t+2]-3*tau[t+1]-4*tau[t-1]-4*tau[t-2]-4*tau[t-3]-4*tau[t-4]-4*tau[t-5]
                                                     -4*tau[t-6]-4*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    t=TIME-3
    tau[t] = ((5/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+4]-2*tau[t+3]-3*tau[t+2]-4*tau[t+1]-5*tau[t-1]-5*tau[t-2]-5*tau[t-3]-5*tau[t-4]-5*tau[t-5]
                                                     -5*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    t = TIME-4
    tau[t] = ((6/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+5]-2*tau[t+4]-3*tau[t+3]-4*tau[t+2]-5*tau[t+1]-6*tau[t-1]-6*tau[t-2]-6*tau[t-3]-6*tau[t-4]-6*tau[t-5]
                                                     -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t])) 
    
    t=TIME-5
    tau[t] = ((7/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+6]-2*tau[t+5]-3*tau[t+4]-4*tau[t+3]-5*tau[t+2]-6*tau[t+1]-7*tau[t-1]-7*tau[t-2]-7*tau[t-3]-7*tau[t-4]-7*tau[t-5]
                                                     -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t])) 
    
    t=TIME-6
    tau[t] = ((8/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+7]-2*tau[t+6]-3*tau[t+5]-4*tau[t+4]-5*tau[t+3]-6*tau[t+2]-7*tau[t+1]-8*tau[t-1]-8*tau[t-2]-8*tau[t-3]-8*tau[t-4]-7*tau[t-5]
                                                     -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t])) 
    
    t=TIME-7
    tau[t] = ((9/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+8]-2*tau[t+7]-3*tau[t+6]-4*tau[t+5]-5*tau[t+4]-6*tau[t+3]-7*tau[t+2]-8*tau[t+1]-9*tau[t-1]-9*tau[t-2]-9*tau[t-3]-8*tau[t-4]-7*tau[t-5]
                                                     -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                          +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t])) 
    
    t=TIME-8
    tau[t] = ((10/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+9]-2*tau[t+8]-3*tau[t+7]-4*tau[t+6]-5*tau[t+5]-6*tau[t+4]-7*tau[t+3]-8*tau[t+2]-9*tau[t+1]-10*tau[t-1]-10*tau[t-2]-9*tau[t-3]-8*tau[t-4]-7*tau[t-5]
                                                      -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))  
    
    t=TIME-9
    tau[t] = ((11/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-tau[t+10]-2*tau[t+9]-3*tau[t+8]-4*tau[t+7]-5*tau[t+6]-6*tau[t+5]-7*tau[t+4]-8*tau[t+3]-9*tau[t+2]-10*tau[t+1]-11*tau[t-1]-10*tau[t-2]-9*tau[t-3]-8*tau[t-4]-7*tau[t-5]
                                                      -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))  
    
    
    for(t in 12:(TIME-10)){
      
      tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-11*tau[t-1]-10*tau[t-2]-9*tau[t-3]-8*tau[t-4]-7*tau[t-5]
                                                        -6*tau[t-6]-5*tau[t-7]-4*tau[t-8]-3*tau[t-9]-2*tau[t-10]-tau[t-11]
                                                        -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]-8*tau[t+4]-7*tau[t+5]
                                                        -6*tau[t+6]-5*tau[t+7]-4*tau[t+8]-3*tau[t+9]-2*tau[t+10]-tau[t+11])
                                             +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    }
    
    
    
    t=11
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-2*tau[t-10]-3*tau[t-9]-4*tau[t-8]-5*tau[t-7]-6*tau[t-6]-7*tau[t-5]-8*tau[t-4]
                                                      -9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=10
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-3*tau[t-9]-4*tau[t-8]-5*tau[t-7]-6*tau[t-6]-7*tau[t-5]-8*tau[t-4]
                                                      -9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=9
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-4*tau[t-8]-5*tau[t-7]-6*tau[t-6]-7*tau[t-5]-8*tau[t-4]
                                                      -9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=8
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-5*tau[t-7]-6*tau[t-6]-7*tau[t-5]-8*tau[t-4]-9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    
    t=7
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-6*tau[t-6]-7*tau[t-5]-8*tau[t-4]-9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=6
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-7*tau[t-5]-8*tau[t-4]-9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    
    t=5
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-8*tau[t-4]-9*tau[t-3]-10*tau[t-2]-11*tau[t-1]
                                                      -11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=4
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-9*tau[t-3]-10*tau[t-2]-11*tau[t-1]-11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=3
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-10*tau[t-2]-11*tau[t-1]-11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=2
    tau[t] = ((12/sig2tau+1/v[t-1])^(-1))*(1/sig2tau*(-11*tau[t-1]-11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                      -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                      -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                      -2*tau[t+10]-tau[t+11])
                                           +1/v[t-1]*(y[t-1]-theta[,t]%*%matrix(X[t-1,])-u[t]))
    
    
    t=1
    tau[t] = ((11/sig2tau+1/C0tau)^(-1))*(m0tau/C0tau+1/sig2tau*(-11*tau[t+1]-10*tau[t+2]-9*tau[t+3]
                                                                 -8*tau[t+4]-7*tau[t+5]-6*tau[t+6]
                                                                 -5*tau[t+7]-4*tau[t+8]-3*tau[t+9]
                                                                 -2*tau[t+10]-tau[t+11]))
    
    
    # STORE SAMPLES
    
    store_gamma[i,,] = t(gamma) 
    store_theta[i,,] = theta 
    store_v[i,] = v
    store_d[i,] = d
    store_u[i,] = u
    store_tau[i,] = tau
    store_reg[i,] = reg
  }
  
  if(missing(new.X)){f = NULL }else{
    
    GG.1 = diag(c(phi1*gamma[dim(gamma)[1],],1,1),P+1+1)
    state_vec = matrix(c(theta[,TIME+1],u[TIME+1],tau[TIME+1]),ncol=1)
    a = GG.1%*%state_vec
    f = FF.1%*%a
    
  }
  
  cat("MCMC takes")
  print(proc.time() - ptm)
  
  cat("")
  
  return(list(beta=store_theta,gamma=store_gamma,
              v=store_v,u=store_u,d=store_d,tau = store_tau,
              fc.f = f,fit.reg=store_reg,N=N,TIME=TIME))
  
}
