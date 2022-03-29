#---------------------------------------------------------------------------#
#                                                                           #
# SSVS with DSS for Dynamic Linear Models (Stochastic Volatility)           #
#                                                                           #
#---------------------------------------------------------------------------#


DSSVS_SV<-function(y,X,N,gamma0,phi0,v0,phi1,THETA,lambda0,lambda1,
                   start_params,start_latent){
  
  cat("Starting BSTS with dynamic Shrinkage")
  
  TIME <- length(y)
  P <- dim(X)[2]
  
  require(dlm)
  require(stochvol)
  
  if(missing(start_params)){
    params <- list(mu = -10, phi = 0.9, sigma = 0.1,
                   nu = Inf, rho = 0, beta = NA,
                   latent0 = -10)
  }else{
    params <-start_params
  }
  
  if(missing(start_latent)){
    start_latent = rep(-10, length(y))
  }else{
    start_latent = start_latent
  }
  
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P,TIME+1))
  store_v <-matrix(NA,N,TIME+1)
  store_f <-matrix(NA,N,TIME)
  
  v<-rep(v0,TIME+1)
  
  #gamma<-matrix(c(rep(1,4),rep(0,16)),TIME+1,P,byrow=T)
  gamma<-matrix(gamma0,TIME+1,P,byrow=T)
  
  ptm <- proc.time()
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  for(i in 1:N){
    
    if(i %% 100==0) cat("\n", i, "loops...") 
    if(i == N) cat("Done!\n") 
    
    X = X
    G = phi1*gamma[2:(TIME+1),]
    W = gamma[2:(TIME+1),]*lambda1+(1-gamma[2:(TIME+1),])*lambda0
    V = v[2:(TIME+1)]
    
    MATX<-cbind(X,G,W,V)
    
    mod<-dlm(m0=phi0*gamma[1,],
             C0=diag(gamma[1,]*lambda1/(1-phi1^2)+(1-gamma[1,])*lambda0),
             FF=matrix(rep(1,P),nrow=1),
             GG=diag(1,P),
             JFF=matrix(c(1:P),nrow=1),
             JGG=diag(c((P+1):(P+P))),
             V=1,
             JV = P+P+P+1,
             W = diag(1,P),
             JW = diag(c((P+P+1):(P+P+P))),
             X=MATX)
    
    Filt<-dlmFilter(y,mod)
    f<-Filt$f
    theta<-t(dlmBSample(Filt))
    
    beta<-matrix(theta,ncol=1,byrow=T)
    thetabeta<-matrix(theta[,2:(TIME+1)],ncol=1,byrow=T)
    
    thetaind = THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                    +(1-THETA)*dnorm(beta[1:((TIME+1)*P-P),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(beta[1:((TIME+1)*P-P),1]-phi0)
    
    pstar = thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))/(thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))
                                                                         +(1-thetaind)*dnorm(beta[(P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
    
    gammino = rbinom(TIME*P,1,pstar)
    
    gammino0 = rbinom(P,1,thetaind[1:P])
    
    gamma = matrix(c(gammino0,gammino),ncol=P,byrow=T)
    
    
    res=y-bigX%*%thetabeta
    
    if(i==1){
      v[1:(TIME+1)]=0.25
    }else{
      
      res_fast = svsample_fast_cpp(res,
                                   startpara = params, startlatent = start_latent)
      
      v[2:(TIME+1)]=exp(res_fast$latent)
      v[1]=exp(res_fast$latent0)
    }
    
    store_gamma[i,,] = t(gamma) 
    store_theta[i,,] = theta
    store_f[i,] = f
    store_v[i,] = v
    
  }
  
  cat("MCMC takes:")
  print(proc.time() - ptm)
  cat("Note: this function used a log-AR(1) process for volatility")
  
  return(list(f=store_f,theta=store_theta,gamma=store_gamma,v=store_v))
  
}
