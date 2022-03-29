# This function performs MAP smoothing in 
# TVP regression models with stochastic volatility using 
# Dynamic EMVS with a Particle Smoothing routine
# -----------------------------------------------



DEMVS_PS<-function(y,X,N,phi0,phi1,THETA,lambda0,lambda1,new.X){
  
  # Particle Filter
  #--------------------------------------------------------
  
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
  #-----------------------------------------------------------------
  
  require(matlib)
  cat("Starting EMVS for dynamic Shrinkage in Dynamic Linear Model")
  
  P<-dim(X)[2]
  TIME<-length(y)
  
  ols = lm(y ~ X)
  beta0 = ols$coefficients[-1]
  beta0[is.na(beta0)]=0
  
  theta<-matrix(beta0,nrow=P,ncol=TIME+1)
  v<-c()
  precision<-c()
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  if(missing(new.X)){}else{
    FF.1 = matrix(c(new.X),nrow=1)
  } 
  
  # Store objects
  #-------------------------------------------------
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P,TIME+1))
  store_v <-matrix(NA,N,TIME)
  
  # Start Gibbs sampling
  #-------------------------------------------------
  
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
    
    if(i<10){
      precision = rep(4,TIME)
      v=1/precision
    }else{
      
      out.pf = BPF_SV(res,N=1000,alpha0=-2,alpha1=0.9,sigzeta=0.1)
      h = out.pf$xhat
      v = exp(h)
      precision = 1/v
    }
    
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
    store_v[i,] = v
    
  }
  
  if(missing(new.X)){f = NULL }else{
    
    GG.1 = diag(phi1*gamma[dim(gamma)[1],],P)
    state_vec = matrix(c(theta[,TIME+1]),ncol=1)
    a = GG.1%*%state_vec
    f = FF.1%*%a
    
  }
  cat("MCMC takes")
  print(proc.time() - ptm)
  
  return(list(beta=store_theta,fc.f = f,gamma=store_gamma,v=store_v))
  
}







bb = BPF_SV(res,N,alpha0,alpha1,sigzeta)

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


x<-rnorm(N,alpha0+alpha1*x+(sigzeta^2)/4*((data[t]^2)*exp(-(alpha0+alpha1*x))-2),sigzeta)
w1<-w*dnorm(data[t],0,exp(x/2))*dnorm(x,alpha0+alpha1*xprev,sigzeta)/dnorm(x,mean=alpha0+alpha1*xprev+(sigzeta^2))
w = w1/sum(w1)



set.seed(12345)
n     =  100
al =  0
be =  0.9
ta  =  0.1
y4     = rep(0,n)
x4     = rep(0,n)
x4[1]  = al/(1-be)
y4[1]  = rnorm(1,0,exp(x4[1]/2))
for (t in 2:n){
  x4[t] = al + be*x4[t-1] + rnorm(1,0,ta)
  y4[t] = rnorm(1,0,exp(x4[t]/2))
}



pippo = particleFilterSVmodel(res,c(-2,0.9,0.1),1000)





svlm = SVLWfun(res,1000,0,100,-2,0.1,0.9,0.1,1,4)
sapply(1:TIME, function(i) weighted.mean(svlm$pars[,1,i],bb$ws[i,]))



applyweighted.mean(svlm$xs[1,],svlm$ws[1,])


SVLWfun<-function(data,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda){
  xs = rnorm(N,m0,sqrt(C0))
  pars = cbind(rnorm(N,ealpha,sqrt(valpha)),rnorm(N,ebeta,sqrt(vbeta)),
               log(1/rgamma(N,nu/2,nu*lambda/2)))
  delta = 0.75
  a = (3*delta-1)/(2*delta)
  h2 = 1-a^2
  parss = array(0,c(N,3,length(data)))
  xss = matrix(NA,nrow=length(data),ncol=N)
  ws = matrix(NA,nrow=length(data),ncol=N)
  w = rep(1/N,N)
  mpar <-c()
  for (t in 1:length(data)){
    for(i in 1:3){
      mpar[i] = weighted.mean(pars[,i],w)
    }
    vpar = var(pars)
    ms = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    mus = pars[,1]+pars[,2]*xs
    weight = w*dnorm(data[t],0,exp(mus/2))
    k = sample(1:N,size=N,replace=T,prob=weight)
    ms1 = ms[k,] + matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w = dnorm(data[t],0,exp(xt/2))/dnorm(data[t],0,exp(mus[k]/2))
    w = w/sum(w)
    ESS = 1/sum(w^2)
    if(ESS<(N/2)){
      index<-sample(N,size=N,replace=T,prob=w)
      xs<-xt[index]
      pars<-ms1[index,]
      w<-rep(1/N,N)
    }else{
      xs<-xt
      pars<-ms1
    }
    xss[t,] = xs
    parss[,,t] = pars
    ws[t,] = w
  }
  for(p in 1:3){
    param[,p] = sapply(1:TIME, function(i) weighted.mean(parss[,p,i],ws[i,]))
  }
  for(t in (length(data)-1):1){
    print(t)
    ws[t,] = dnorm(xss[t,],param[t,1]+param[t,2]*(xss[t+1,]-alpha0),param[t,3])
    index<-sample(N,size=N,replace=T,prob=ws[t,])
    xss[t,]<-xss[t,index]
  }
  return(list(xs=xss,pars=parss,ws=ws,vpars=vpars))
}


