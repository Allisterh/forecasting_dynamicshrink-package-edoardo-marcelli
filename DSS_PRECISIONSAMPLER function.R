# Dynamic SSVS (precision sampler) for TVP regression models with SV
#--------------------------------------------------------------------

DSS_PRECISIONSAMPLER = function(y,X,N,burn,THETA,lambda1,
                                lambda0,phi1,phi0,gamma0,v0,
                                start_params,start_latent,new.X,activate=F){
  
  # Packages
  
  require(Matrix)
  require(stochvol)
  require(pracma)
  require(Rcpp)
  require(MASS)
  require(matrixStats)
  require(resample)
  cppFunction(plugins=c("cpp11"), "NumericVector cpprbinom2(int n, double size, NumericVector prob) { 
    NumericVector v = no_init(n);
    std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); }); 
    return(v);}")
  
  # Stocvol Package options
  
  if(missing(start_params)){
    params <- list(mu = -1000, phi = 0.9, sigma = 0.1,
                   nu = Inf, rho = 0, beta = NA,
                   latent0 = -1000)
  }else{
    params <-start_params
  }
  
  if(missing(start_latent)){
    start_latent = rep(-1000, length(y))
  }else{
    start_latent = start_latent
  }
  
  # Preliminaries
  
  TIME <- length(y)
  P <- dim(X)[2]
  bigX <- Matrix(nrow = TIME, ncol = TIME*P, data = 0, sparse = TRUE)
  j=1
  for(t in 1:TIME){
    bigX[t,j:(j+P-1)]<-X[t,]
    j=j+P
  }
  
  store_gammino = matrix(NA,N,TIME*P)
  store_gammino0 = matrix(NA,N,P)
  store_Sigtheta = matrix(NA,N,TIME*P)
  store_Sigtheta0 = matrix(N,1)
  store_theta = matrix(NA,N,TIME*P)
  store_theta0 = matrix(NA,N,P)
  store_h = matrix(NA,N,TIME)
  step_ahead=1
  store_fc.m<-matrix(NA,nrow=step_ahead,ncol=N)
  store_fc.y<-matrix(NA,nrow=step_ahead,ncol=N)
  store_fc.v<-matrix(NA,nrow=step_ahead,ncol=N)
  
  h = matrix(log(v0),TIME,1)
  theta0 = matrix(0,P,1)
  gammino = gamma0*matrix(1,TIME*P,1)
  gammino0 = gamma0*matrix(1,P,1)
  
  if(missing(new.X)){}else{
    FF.new = matrix(new.X,nrow=1)
  }
  
  Htheta = Matrix(nrow = TIME*P, ncol = TIME*P, data = 0, sparse = TRUE)
  diag(Htheta)=1
  Htheta = Htheta - sparseMatrix(i=c((P+1):(TIME*P)),j=c(1:((TIME-1)*P)),
                                 x=c(phi1*gammino[(P+1):length(gammino)]),dims=c(TIME*P,TIME*P))
  
  Sigtheta = sparseMatrix(i=c(1:(TIME*P)),j=c(1:(TIME*P)),
                          x=c(gammino*lambda1+(1-gammino)*lambda0),dims=c(TIME*P,TIME*P))
  
  Sigtheta0 = gammino0*lambda1+(1-gammino0)*lambda0
  
  invSig = Diagonal(TIME)
  invS = Diagonal(TIME*P)
  a0 = Matrix(nrow=TIME*P,ncol=1,data=0,sparse=T)
  Ktheta0 = Diagonal(P)
  I_tp=Diagonal(TIME*P)
  
  # START GIBBS SAPLING
  
  ptm <- proc.time()
  
  for(it in 1:(N+burn)){
    
    if(it %% 100==0) cat("\n", it, "loops...") 
    if(it == N) cat("Done!\n") 
    
    # Sample states 1:T
    
    diag(invSig)=c(exp(-h))
    
    vec = 1/diag(Sigtheta)
    diag(invS) = vec
    
    XinvSig = t(bigX)%*%invSig
    
    HinvSH = t(Htheta)%*%invS%*%Htheta
    
    a0[1:P] =  phi1*gammino[1:P,1]*theta0
    
    alptheta = solve(Htheta,a0)
    
    Ktheta = HinvSH + XinvSig%*%bigX
    
    dtheta = XinvSig%*%y + HinvSH%*%alptheta
    
    thetahat = solve(Ktheta,dtheta)
    
    theta = thetahat + solve(t(chol(Ktheta)),rnorm(TIME*P))
    
    # Sample state 0
    
    vec = 1/diag(Sigtheta0)
    diag(Ktheta0) = vec
    
    a0.0 = theta[1:P,]/Sigtheta0
    theta0hat = solve(Ktheta0,a0.0)
    
    theta0 = theta0hat + solve(t(chol(Ktheta0)),rnorm(P))
    
    # Sample indicators
    
    
    thetaind = THETA*dnorm(theta[(P+1-P):(TIME*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(theta[(P+1-P):(TIME*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))+
                                                                                         (1-THETA)*dnorm(theta[(P+1-P):(TIME*P-P),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(theta[(P+1-P):(TIME*P-P),1]-phi0)
    
    
    pstar = thetaind*dnorm(theta[(P+1):(TIME*P),1],mu,sqrt(lambda1))/(thetaind*dnorm(theta[(P+1):(TIME*P),1],mu,sqrt(lambda1))+
                                                                        (1-thetaind)*dnorm(theta[(P+1):(TIME*P),1],0,sqrt(lambda0)))
    
    pstar[is.na(pstar)]=1
    gammino[(P+1):(TIME*P),1] = cpprbinom2(TIME*P-P,1,pstar)
    
    
    thetaind = THETA*dnorm(theta0[1:P,1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(theta0[1:P,1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                         +(1-THETA)*dnorm(theta0[1:P,1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(theta0[1:P,1]-phi0)
    
    pstar = thetaind*dnorm(theta[1:P,1],mu,sqrt(lambda1))/(thetaind*dnorm(theta[1:P,1],mu,sqrt(lambda1))
                                                           +(1-thetaind)*dnorm(theta[1:P,1],0,sqrt(lambda0)))
    
    pstar[is.na(pstar)]=1
    gammino[1:P,1] = cpprbinom2(P,1,pstar)
    thetaind[is.na(thetaind)]=1
    
    gammino0[1:P,1] = cpprbinom2(P,1,thetaind)
    
    if(activate==T){
      gammino[seq(1,TIME*P,by=P),] = 1
    }
    
    # Compute state variance matrix
    
    diag(Sigtheta) = c(gammino*lambda1+(1-gammino)*lambda0)
    
    Sigtheta0 = gammino0*lambda1+(1-gammino0)*lambda0
    
    # Compute H theta
    
    Htheta = I_tp - sparseMatrix(i=c((P+1):(TIME*P)),j=c(1:((TIME-1)*P)),
                                 x=c(phi1*gammino[(P+1):length(gammino)]),dims=c(TIME*P,TIME*P))
    
    # Stocvol
    
    res=y-bigX%*%theta
    
    if(it<10){
      h=rep(-10,TIME)
    }else{
      res_fast = svsample_fast_cpp(as.matrix(res),
                                   startpara = params, startlatent = start_latent)
      
      h=res_fast$latent
      alpha0 = res_fast$para[1,1]
      alpha1 = res_fast$para[1,2]
      sigzeta = res_fast$para[1,3]
      h.last = res_fast$latent[1,TIME]
      h.new = alpha0 + alpha1*h.last + sigzeta*rnorm(1,0,1)
    }
    
    
    if(it > burn){
      store_theta[it-burn,] = t(theta)[,1:(TIME*P)]
      store_h[it-burn,] = t(h)
      store_theta0[it-burn,] = t(theta0)[1:P]
      store_gammino[it-burn,] = t(gammino)
      store_gammino0[it-burn,] = t(gammino0)
      
      if(missing(new.X)){}else{
        thetaind.new = THETA*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                                        +(1-THETA)*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],0,sqrt(lambda0)))
        thetaind[is.na(thetaind.new)]=1
        gammino.new = cpprbinom2(P,1,thetaind.new)
        
        WW.new = diag(gammino.new*lambda1+(1-gammino.new)*lambda0)
        V.new = exp(h.new)
        
        m.new = thetahat[(TIME*P-P+1):(TIME*P),1]
        GG.new = diag(phi1*gammino.new)
        a.new = GG.new%*%m.new
        f.new = FF.new%*%a.new
        
        C.new = solve(Ktheta[(TIME*P-P+1):(TIME*P),(TIME*P-P+1):(TIME*P)])
        R.new = (GG.new)%*%C.new%*%t(GG.new)+WW.new
        Q.new = FF.new%*%R.new%*%t(FF.new)+V.new
        Q.new.mat = as.matrix(Q.new)
        
        store_fc.m[,it-burn]<-f.new
        store_fc.v[,it-burn]<-Q.new.mat
      }
      
    }else{}
    
  }
  
  cat("MCMC takes:")
  print(proc.time() - ptm)
  
  est.theta = colMeans(store_theta)
  est.h = colMedians(store_h)
  est.theta0 = colMeans(store_theta0)
  est.gammino = colMeans(store_gammino)
  est.gammino0 = colMeans(store_gammino0)
  
  est.var.theta = colVars(store_theta)
  est.var.theta0 = colVars(store_theta0)
  
  theta = matrix(est.theta,nrow=TIME,ncol=P,byrow=T)
  h = matrix(est.h,nrow=TIME,ncol=1,byrow=T)
  theta0 = matrix(est.theta0,nrow=1,ncol=P,byrow=T)
  gammino = matrix(est.gammino,nrow=TIME,ncol=P,byrow=T)
  gammino0 = matrix(est.gammino0,nrow=TIME,ncol=P,byrow=T)
  var.theta = matrix(est.var.theta,nrow=TIME,ncol=P,byrow=T)
  var.theta0 = matrix(est.var.theta0,nrow=1,ncol=P,byrow=T)
  
  est.fcm = median(store_fc.m[1,])
  est.fcv = median(store_fc.v[1,])
  
  
  return(list(beta = theta, h=h,beta0 = theta0, ind=gammino, ind0=gammino0,
              beta.sample=store_theta, h.sample=store_h,
              var.beta=var.theta,var.beta0=var.theta0,
              fc.f = est.fcm, fc.Q = est.fcv
              ,beta0.sample=store_theta0,
              gammino.sample=store_gammino,gammino0.sample=store_gammino,
              bigX=bigX,y=y,
              THETA=THETA,lambda0=lambda0,lambda1=lambda1,TIME=TIME,X=X,phi0=phi0,
              phi1=phi1))
} 