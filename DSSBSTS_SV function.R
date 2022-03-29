# Dynamic shrinkage in BSTS models and SV 
#------------------------------------------

DSSBSTS_SV<-function(y,X,N,S,U,gamma0,phi0,v0,phi1,THETA,lambda0,lambda1,
                     sig2u=0.1,sig2d=0.1,sig2tau=0.1,
                     start_params,start_latent,new.X,activate=F){
  
  cat("Starting BSTS with dynamic Shrinkage")
  
  # Preliminaries
  #--------------
  
  TIME <- length(y)
  P <- dim(X)[2]
  z<-X
  
  require(dlm)
  require(stochvol)
  require(Rcpp)
  require(MASS)
  cppFunction(plugins=c("cpp11"), "NumericVector cpprbinom2(int n, double size, NumericVector prob) { 
    NumericVector v = no_init(n);
    std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); }); 
    return(v);}")
  
  if(missing(start_params)){
    params <- list(mu = -2, phi = 0.9, sigma = 0.1,
                   nu = Inf, rho = 0, beta = NA,
                   latent0 = -2)
  }else{
    params <-start_params
  }
  
  if(missing(start_latent)){
    start_latent = rep(-2, length(y))
  }else{
    start_latent = start_latent
  }
  
  if(missing(S)){S=0}else{
    S=S
    Seas<-matrix(c(1,rep(0,S-1)),nrow=TIME,ncol=S,byrow=T)}
  
  if(missing(U)){U=0}else{
    U=U
    Lev<-matrix(c(1,0),nrow=TIME,ncol=U,byrow=T)}
  
  # Kalman recursions matrix
  #-------------------------
  FF.1 = matrix(rep(1,P),nrow=1)
  GG.1 = diag(1,P)
  JFF.1 = matrix(c(1:P),nrow=1)
  JGG.1 = diag(c((P+1):(P+P)))
  V.1 = 1
  JV.1 = P+P+P+1
  W.1 = diag(1,P)
  JW.1 = diag(c((P+P+1):(P+P+P)))
  
  # Store objects
  #---------------
  
  store_gamma<-array(NA,c(N,P,TIME+1))
  store_theta <-array(NA,c(N,P+S+U,TIME+1))
  store_v <-matrix(NA,N,TIME)
  store_reg <-matrix(NA,N,TIME)
  store_y.new <-matrix(NA,N,1)
  store_V.new <-matrix(NA,N,1)
  store_q <-matrix(NA,N,TIME)
  store_c <-array(NA,c(N,P+S+U,TIME))
  
  # Start volatility and indicators
  #--------------------------------
  
  v<-rep(v0,TIME)
  gamma<-matrix(gamma0,TIME+1,P,byrow=T)
  
  # Prepare other usefull stuff
  #----------------------------
  
  if(S==0 & U!=0){
    augX<-cbind(X,Lev)
  }else if(S!=0 & U==0){
    augX<-cbind(X,Seas)
  }else if(S==0 & U==0){
    augX<-X
  }else{
    augX<-cbind(X,Seas,Lev)
  }
  
  if(missing(new.X)){}else{
    
    if(S==0 & U!=0){
      F_1<-c(new.X,Lev[1,])
    }else if(S!=0 & U==0){
      F_1<-c(new.X,Seas[1,])
    }else if(S==0 & U==0){
      F_1<-new.X
    }else{
      F_1<-c(new.X,Seas[1,],Lev[1,])
    }
    
  }
  
  bigX<-matrix(0,nrow=TIME,ncol=TIME*(P+S+U))
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P+S+U-1)]<-augX[t,]
    j=j+P+S+U
  }
  
  smallX<-matrix(0,nrow=TIME,ncol=TIME*P)
  j=1
  for(t in 1:TIME){
    smallX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  # Start Gibbs sampling
  #---------------------
  
  ptm <- proc.time()
  
  for(i in 1:N){
    
    if(i %% 100==0) cat("\n", i, "loops...") 
    if(i == N) cat("Done!\n") 
    
    X = X
    G = phi1*gamma[2:(TIME+1),]
    W = gamma[2:(TIME+1),]*lambda1+(1-gamma[2:(TIME+1),])*lambda0
    V = matrix(v)
    
    MATX<-cbind(X,G,W,V)
    
    mod<-dlm(m0=phi0*gamma[1,],
             C0=diag(gamma[1,]*lambda1/(1-phi1^2)+(1-gamma[1,])*lambda0),
             FF=FF.1,
             GG=GG.1,
             JFF=JFF.1,
             JGG=JGG.1,
             V=V.1,
             JV = JV.1,
             W = W.1,
             JW = JW.1,
             X=MATX)
    
    if(S==0 & U!=0){
      mod = mod+dlmModPoly(order = U,dV=0,dW = c(sig2u, sig2d))
    }else if(S!=0 & U==0){
      mod = mod+dlmModSeas(S+1,dV=0,dW = c(sig2tau, rep(0, S+1-2)))
    }else if(S==0 & U==0){
      mod = mod
    }else{
      mod = mod+dlmModSeas(S+1,dV=0,dW = c(sig2tau, rep(0, S+1-2)))+dlmModPoly(order = U,dV=0,dW = c(sig2u, sig2d))
    }
    
    Filt<-dlmFilter(y,mod)
    
    theta<-t(dlmBSample(Filt))
    
    beta<-matrix(theta[1:P,],ncol=1,byrow=T)
    thetabeta<-matrix(theta[,2:(TIME+1)],ncol=1,byrow=T)
    
    thetaind = THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                    +(1-THETA)*dnorm(beta[1:((TIME+1)*P-P),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(beta[1:((TIME+1)*P-P),1]-phi0)
    
    pstar = thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))/(thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))
                                                                         +(1-thetaind)*dnorm(beta[(P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
    
    gammino = cpprbinom2(TIME*P,1,pstar)
    
    gammino0 = cpprbinom2(P,1,thetaind[1:P])
    
    gamma = matrix(c(gammino0,gammino),ncol=P,byrow=T)
    
    if(activate==T){
      gamma[,1]=1 
    }else{}
    
    res=y-bigX%*%thetabeta
    reg = smallX%*%beta[(P+1):((TIME+1)*P),1]
    
    if(i<10){
      v=v
      h.new = log(v0)
    }else{
      res_fast = svsample_fast_cpp(res, burnin=100,
                                   startpara = params, startlatent = start_latent)
      
      v=exp(res_fast$latent)
      alpha0 = res_fast$para[1,1]
      alpha1 = res_fast$para[1,2]
      sigzeta = res_fast$para[1,3]
      h.last = res_fast$latent[1,TIME]
      h.new = alpha0 + alpha1*h.last + sigzeta*rnorm(1,0,1)
    }
    
    store_gamma[i,,] = t(gamma) 
    store_theta[i,,] = theta
    store_reg[i,] = reg
    store_v[i,] = v
    
    if(missing(new.X)){}else{
      
      # KALMAN FORECAST
      thetaind.new = THETA*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                                         +(1-THETA)*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
      
      gammino.new = cpprbinom2(P,1,thetaind.new)
      
      V.new = rbind(V,exp(h.new))
      G.new = rbind(G,phi1*gammino.new)
      W.new = rbind(W,
                    gammino.new*lambda1+(1-gammino.new)*lambda0)
      X_new = rbind(X,new.X)
      
      MATX<-cbind(X_new,G.new,W.new,V.new)
      
      mod<-dlm(m0=phi0*gamma[1,],
               C0=diag(gamma[1,]*lambda1/(1-phi1^2)+(1-gamma[1,])*lambda0),
               FF=FF.1,
               GG=GG.1,
               JFF=JFF.1,
               JGG=JGG.1,
               V=V.1,
               JV = JV.1,
               W = W.1,
               JW = JW.1,
               X=MATX)
      
      if(S==0 & U!=0){
        mod = mod+dlmModPoly(order = U,dV=0,dW = c(sig2u, sig2d))
      }else if(S!=0 & U==0){
        mod = mod+dlmModSeas(S+1,dV=0,dW = c(sig2tau, rep(0, S+1-2)))
      }else if(S==0 & U==0){
        mod = mod
      }else{
        mod = mod+dlmModSeas(S+1,dV=0,dW = c(sig2tau, rep(0, S+1-2)))+dlmModPoly(order = U,dV=0,dW = c(sig2u, sig2d))
      }
      
      y1 = c(y,y[length(y)])
      
      Forc<-dlmFilter(y1,mod)
      listR <- dlmSvd2var(Forc$U.R, Forc$D.R)
      R_DLM = listR[[TIME+1]]
      Q_DLM = F_1%*%R_DLM%*%t(t(F_1))+exp(h.new)
      
      store_y.new[i,]<-Forc$f[length(Forc$f)]
      store_V.new[i,]<-Q_DLM
    }
    
  }
  
  cat("MCMC takes:")
  print(proc.time() - ptm)
  
  return(list(fit.reg=store_reg,theta=store_theta,gamma=store_gamma,v=store_v,bigX=bigX,y=y,
              THETA=THETA,lambda0=lambda0,lambda1=lambda1,S=S,U=U,TIME=TIME,X=X,phi0=phi0,
              phi1=phi1,fc.f=store_y.new, fc.Q = store_V.new))
  
}
