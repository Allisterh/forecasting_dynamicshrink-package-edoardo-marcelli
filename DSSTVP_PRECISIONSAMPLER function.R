# Dynamic SSVS (Precision Sampler) for TVP-VAR models with SV
#-------------------------------------------------------------

DSSTVP_PRECISIONSAMPLER = function(X,constant,lags,N,burn,THETA,lambda1,lambda0,phi1,phi0,gamma0,v0,start_params,start_latent){
  
  # Packages
  
  require(Matrix)
  require(stochvol)
  require(pracma)
  require(Rcpp)
  require(MASS)
  require(matlib)
  require(mvtnorm)
  require(matrixStats)
  cppFunction(plugins=c("cpp11"), "NumericVector cpprbinom2(int n, double size, NumericVector prob) { 
    NumericVector v = no_init(n);
    std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); }); 
    return(v);}")
  require(bvarsv)
  list <- unclass(lsf.str(envir = asNamespace("bvarsv"), all = T))[]
  for(i in list){
    assign(paste0(i),eval(parse(text = paste('bvarsv:::', i, sep=''))))
  }
  
  
  # Preliminaries
  
  TIME = dim(X)[1]
  n = dim(X)[2]
  h=lags
  m = n*(n-1)/2
  P = n^2*h
  step_ahead=1
  np=sum(1:n)
  
  XX<-X[(h+1):TIME,]
  Y = matrix(t(XX),ncol=1)
  bigY<-matrix(NA,nrow=TIME-h,ncol=n*h)
  j=1
  for(i in h:1){
    bigY[,j:(j+n-1)]<-X[i:(TIME-h+i-1),]
    j=j+n  
  }
  
  bigY.new<-matrix(NA,nrow=TIME-h+1,ncol=n*h)
  j=1
  for(i in h:1){
    bigY.new[,j:(j+n-1)]<-X[i:(TIME-h+i),]
    j=j+n  
  }
  
  TIME = dim(bigY)[1]
  
  if(constant==TRUE){
    bigY = cbind(matrix(1,nrow=dim(bigY)[1],ncol=1),bigY)
    bigY.new = cbind(matrix(1,nrow=dim(bigY.new)[1],ncol=1),bigY.new)
    pars = (n^2)*h+n
    P = pars
  }else{}
  
  ff.new = bigY.new[dim(bigY.new)[1],]
  FF.new = kronecker(diag(1,n),matrix(ff.new,nrow=1))
  
  Fill_bigX<-list()
  for(i in 1:dim(bigY)[1]){
    block_bigX = t(kronecker(diag(1,n),bigY[i,]))
    Fill_bigX[[i]] = block_bigX 
  }
  bigX = .bdiag(Fill_bigX)
  
  # Stocvol Package options
  
  if(missing(start_params)){
    params <- list(mu = -2, phi = 0.9, sigma = 0.1,
                   nu = Inf, rho = 0, beta = NA,
                   latent0 = -2)
    
  }else{
    params <-start_params
  }
  
  if(missing(start_latent)){
    start_latent = rep(-2, TIME)
  }else{
    start_latent = start_latent
  }
  
  #Store
  
  store_gammino = matrix(NA,N,TIME*P)
  store_gammino0 = matrix(NA,N,P)
  store_Sigtheta = matrix(NA,N,TIME*P)
  store_Sigtheta0 = matrix(N,1)
  store_theta = matrix(NA,N,TIME*P)
  store_theta0 = matrix(NA,N,P)
  store_h = matrix(NA,N,TIME*n)
  store_fc.m<-array(NA,c(n,step_ahead,N))
  store_fc.y<-array(NA,c(n,step_ahead,N))
  store_fc.v<-array(NA,c(np,step_ahead,N))
  
  # Initialize
  
  h = matrix(log(v0),TIME*n,1)
  h1 = matrix(h,ncol=n,byrow=T)
  h.new = matrix(NA,nrow=1,ncol=n)
  theta0 = matrix(0,P,1)
  gammino = gamma0*matrix(1,TIME*P,1)
  gammino0 = gamma0*matrix(1,P,1)
  
  Htheta = Matrix(nrow = TIME*P, ncol = TIME*P, data = 0, sparse = TRUE)
  diag(Htheta)=1
  Htheta = Htheta - sparseMatrix(i=c((P+1):(TIME*P)),j=c(1:((TIME-1)*P)),
                                 x=c(phi1*gammino[(P+1):length(gammino)]),dims=c(TIME*P,TIME*P))
  
  Sigtheta = sparseMatrix(i=c(1:(TIME*P)),j=c(1:(TIME*P)),
                          x=c(gammino*lambda1+(1-gammino)*lambda0),dims=c(TIME*P,TIME*P))
  
  Sigtheta0 = gammino0*lambda1+(1-gammino0)*lambda0
  
  invSig = Diagonal(TIME*n)
  invS = Diagonal(TIME*P)
  a0 = Matrix(nrow=TIME*P,ncol=1,data=0,sparse=T)
  Ktheta0 = Diagonal(P)
  I_tp=Diagonal(TIME*P)
  
  # START GIBBS SAPLING
  
  ptm <- proc.time()
  
  for(it in 1:(N+burn)){
    
    if(it %% 100==0) cat("\n", it, "loops...") 
    if(it == N+burn) cat("Done!\n") 
    
    # Sample states 1:T
    
    diag(invSig)=c(exp(-h))
    
    vec = 1/diag(Sigtheta)
    diag(invS) = vec
    
    XinvSig = t(bigX)%*%invSig
    
    HinvSH = t(Htheta)%*%invS%*%Htheta
    
    a0[1:P] =  phi1*gammino[1:P,1]*theta0
    
    alptheta = solve(Htheta,a0)
    
    Ktheta = HinvSH + XinvSig%*%bigX
    
    dtheta = XinvSig%*%Y + HinvSH%*%alptheta
    
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
    
    pstar[is.na(pstar)]=0
    gammino[(P+1):(TIME*P),1] = cpprbinom2(TIME*P-P,1,pstar)
    
    
    thetaind = THETA*dnorm(theta0[1:P,1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(theta0[1:P,1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                         +(1-THETA)*dnorm(theta0[1:P,1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(theta0[1:P,1]-phi0)
    
    pstar = thetaind*dnorm(theta[1:P,1],mu,sqrt(lambda1))/
      (thetaind*dnorm(theta[1:P,1],mu,sqrt(lambda1))+(1-thetaind)*dnorm(theta[1:P,1],0,sqrt(lambda0)))
    
    pstar[is.na(pstar)]=0
    
    gammino[1:P,1] = cpprbinom2(P,1,pstar)
    
    thetaind[is.na(thetaind)]=0
    
    gammino0[1:P,1] = cpprbinom2(P,1,thetaind)
    
    # Compute state variance matrix
    
    diag(Sigtheta) = c(gammino*lambda1+(1-gammino)*lambda0)
    
    Sigtheta0 = gammino0*lambda1+(1-gammino0)*lambda0
    
    # Compute H theta
    
    Htheta = I_tp - sparseMatrix(i=c((P+1):(TIME*P)),j=c(1:((TIME-1)*P)),
                                 x=c(phi1*gammino[(P+1):length(gammino)]),dims=c(TIME*P,TIME*P))
    
    # Stocvol
    
    res = Y-bigX%*%theta
    res1 = matrix(res,ncol=n,byrow=T)
    
    if(it<10){
      h=rep(-10,TIME*n)
      h.new = matrix(-10,nrow=1,ncol=n)
    }else{
      
      for(i in 1:n){
        res_fast = svsample_fast_cpp(res1[,i],
                                     startpara = params, startlatent = start_latent)
        
        h1[,i] = res_fast$latent
        sigzeta = res_fast$para[1,3]
        alpha1h = res_fast$para[1,2]
        alpha0h = res_fast$para[1,1]
        h.new[,i] = alpha0h+alpha1h*h1[TIME,i]+rnorm(1,0,sigzeta)
      }
      h = matrix(t(h1),ncol=1)
    }
    
    if(it > burn){
      store_theta[it-burn,] = t(theta)[,1:(TIME*P)]
      store_h[it-burn,] = t(h)
      store_theta0[it-burn,] = t(theta0)[1:P]
      store_gammino[it-burn,] = t(gammino)
      store_gammino0[it-burn,] = t(gammino0)
      
      # FORECAST
      
      thetaind.new = THETA*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                                      +(1-THETA)*dnorm(theta[(((TIME)*P)-P+1):((TIME)*P),1],0,sqrt(lambda0)))
      
      gammino.new = cpprbinom2(P,1,thetaind.new)
      
      WW.new = diag(gammino.new*lambda1+(1-gammino.new)*lambda0)
      V.new = diag(c(exp(h.new)))
      
      m.new = thetahat[(TIME*P-P+1):(TIME*P),1]
      GG.new = diag(phi1*gammino.new)
      a.new = GG.new%*%m.new
      f.new = FF.new%*%a.new
      
      C.new = solve(Ktheta[(TIME*P-P+1):(TIME*P),(TIME*P-P+1):(TIME*P)])
      R.new = (GG.new)%*%C.new%*%t(GG.new)+WW.new
      Q.new = FF.new%*%R.new%*%t(FF.new)+V.new
      Q.new.mat = as.matrix(Q.new)
      
      store_fc.m[,,it-burn]<-f.new
      store_fc.v[,,it-burn]<-vechC(Q.new.mat)
      store_fc.y[,,it-burn]<-t(rmvnorm(1,mean=f.new,sigma= Q.new.mat))
      
      
    }else{}
    
  }
  
  cat("MCMC takes:")
  print(proc.time() - ptm)
  
  est.theta = colMeans(store_theta)
  est.h = colMedians(store_h)
  est.theta0 = colMedians(store_theta0)
  est.gammino = colMeans(store_gammino)
  est.gammino0 = colMeans(store_gammino0)
  
  
  
  est.fcy = rowMedians(store_fc.y[,1,])
  est.fcm = rowMedians(store_fc.m[,1,])
  est.fcv = rowMedians(store_fc.v[,1,])
  
  theta = matrix(est.theta,nrow=TIME,ncol=P,byrow=T)
  h = matrix(est.h,nrow=TIME,ncol=n,byrow=T)
  theta0 = matrix(est.theta0,nrow=TIME,ncol=P,byrow=T)
  gammino = matrix(est.gammino,nrow=TIME,ncol=P,byrow=T)
  gammino0 = matrix(est.gammino0,nrow=TIME,ncol=P,byrow=T)
  
  
  return(list(beta = theta/10, h=h,beta0 = theta0/10, ind=gammino, ind0=gammino0,
              beta.sample=store_theta, h.sample=store_h
              ,beta0.sample=store_theta0,
              gammino.sample=store_gammino,gammino0.sample=store_gammino,
              bigX=bigX,Y=Y,
              THETA=THETA,lambda0=lambda0,lambda1=lambda1,TIME=TIME,X=X,phi0=phi0,
              phi1=phi1,
              fc.m = est.fcm/10, fc.y = est.fcy/10, fc.v = est.fcv/100
  ))
}  