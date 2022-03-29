# recommended function
#--------------------------

# Dynamic Shrinkage in TVP-SVAR model with SV
#---------------------------------------------

DSSTVPVAR_SV<-function(X,constant,N,lags,step_ahead,gamma0,phi0,v0,cov0,phi1,THETA.A,THETA.B,lambda0,lambda1,
                       start_params,start_latent){
  
  cat("Starting: TVP-VAR (Primiceri, 2005) with dynamic Spike-and-Slab process priors")
  
  # PRELIMINARIES 
  # --------------------------------------------------------------------------
  
  matsplitter<-function(M, r, c) {
    rg <- (row(M)-1)%/%r+1
    cg <- (col(M)-1)%/%c+1
    rci <- (rg-1)*max(cg) + cg
    N <- prod(dim(M))/r/c
    cv <- unlist(lapply(1:N, function(x) M[rci==x]))
    dim(cv)<-c(r,c,N)
    cv
  } 
  require(MASS)
  require(mvtnorm)
  bdiag_m <- function(lmat) {
    ## Copyright (C) 2016 Martin Maechler, ETH Zurich
    if(!length(lmat)) return(new("dgCMatrix"))
    stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
              (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
              all(vapply(lmat, dim, integer(2)) == k)) # all of them
    N <- length(lmat)
    if(N * k > .Machine$integer.max)
      stop("resulting matrix too large; would be  M x M, with M=", N*k)
    M <- as.integer(N * k)
    ## result: an   M x M  matrix
    new("dgCMatrix", Dim = c(M,M),
        ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
        i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
        p = k * 0L:M,
        x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
  }
  
  require(Rcpp)
  cppFunction(plugins=c("cpp11"), "NumericVector cpprbinom2(int n, double size, NumericVector prob) { 
    NumericVector v = no_init(n);
    std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); }); 
    return(v);}")
  
  reshape.covariance<-function(array){
    takelowertriangular <- function(x){
      x.out = x[lower.tri(x)==T]
      x.out
    }
    mat = apply(array,3,takelowertriangular)
    tmat = t(mat)
    tmat
  }
  
  require(dlm)
  require(stochvol)
  require(Matrix)
  require(matlib)
  require(bvarsv)
  
  list <- unclass(lsf.str(envir = asNamespace("bvarsv"), all = T))[]
  for(i in list){
    assign(paste0(i),eval(parse(text = paste('bvarsv:::', i, sep=''))))
  }
  
  
  # PREPARATION
  #----------------------------------------------------------------------
  
  n = ncol(X)
  h = lags
  
  TIME = dim(X)[1]
  
  Y<-X[(h+1):TIME,]*10
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
  
  if(constant==T){
    ols = lm(Y ~ bigY)
    m0 = c(ols$coefficients)}else{
      ols = lm(Y ~ bigY)
      m0 = c(ols$coefficients[-1,])
    }
  
  # Stoch vol options
  #-------------------------------------------------------------
  if(missing(start_params)){
    
    params <- list(mu = -2, phi = 0.8, sigma = 0.1,
                   nu = Inf, rho = 0, beta = NA,
                   latent0 = -2)
  }else{
    params <-start_params
  }
  
  if(missing(start_latent)){
    start_latent = rep(-2, dim(Y)[1])
  }else{
    start_latent = start_latent
  }
  
  # COnstant options
  #---------------------------------------------------------------
  if(constant==TRUE){
    bigY = cbind(matrix(1,nrow=dim(bigY)[1],ncol=1),bigY)
    bigY.new = cbind(matrix(1,nrow=dim(bigY.new)[1],ncol=1),bigY.new)
    nh = n*h+1
    pars = (n^2)*h+n
    P = pars
  }else{
    nh = n*h
    pars = (n^2)*h
    P = pars
  }
  
  ff.new = bigY.new[dim(bigY.new)[1],]
  FF.new = kronecker(diag(1,n),matrix(ff.new,nrow=1))
  
  PA = sum(1:(n-1))
  
  TIME = dim(Y)[1]
  
  variance = matrix(v0,TIME,ncol=n)
  covariance =  matrix(cov0,TIME,ncol=PA)
  gamma = matrix(gamma0,TIME+1,P,byrow=T)
  bigSigma = matrix(v0,nrow=TIME,ncol=n)
  h.new = matrix(-2,nrow=1,ncol=n)
  sigtdraw = matrix(0.1,nrow=n,ncol=TIME)
  A.new = diag(1,n,n)
  
  oldY = matrix(t(Y),ncol=1,nrow=n*TIME)
  
  Fill_bigX<-list()
  for(i in 1:dim(bigY)[1]){
    block_bigX = t(kronecker(diag(1,n),bigY[i,]))
    Fill_bigX[[i]] = block_bigX 
  }
  bigX = as.matrix(bdiag(Fill_bigX))
  
  # MATRIX OF FFBS (1)
  #----------------------------------------------------------
  
  FF.1 = kronecker(diag(1,n),matrix(rep(1,nh),nrow=1))
  JFF.1 = kronecker(diag(1,n),matrix(seq(1,nh),nrow=1))
  GG.1 = diag(1,pars)
  JGG.1 = diag(c((nh+1):(nh+pars)))
  W.1 = diag(1,pars)
  JW.1 =diag(c((nh+pars+1):(nh+pars+pars)))
  V.1 = diag(1,n)
  JV.1 = diag(c((nh+pars+pars+1):(nh+pars+pars+n)))
  JV.1[lower.tri(JV.1,diag=F)] <- (nh+pars+pars+n+1):(nh+pars+pars+n+PA)
  
  # MATRIX OF FFBS (2)
  #---------------------------------------------------------------------
  
  FF.2 = matrix(0, nrow=n, ncol=PA)
  k=0
  j=1
  for(i in 1:(n-1)){
    FF.2[2+k,j:(j+length(seq(1:i))-1)] = rep(1,i)
    k = k+1
    j = j+length(seq(1:i))
  }
  
  JFF.2 = matrix(0, nrow=n, ncol=PA)
  k=0
  j=1
  for(i in 1:(n-1)){
    JFF.2[2+k,j:(j+length(seq(1:i))-1)] = seq(1:i)
    k=k+1
    j=j+length(seq(1:i))
  }
  np = sum(FF.2==1)
  GG.2 = diag(1,PA)
  JGG.2 = diag(c(n:(n+np-1)),np)
  W.2 = diag(1,PA)
  JW.2 = diag(c((n+np):(n+2*np-1)),np)
  V.2 = diag(1,n)
  JV.2 = diag(c((n+2*np):(n+2*np+n-1)))
  
  aux = matrix(gamma0,TIME+1,np,byrow=T)
  
  # STORING OBJECTS
  #------------------------------------------------------------
  
  store_gamma<-array(NA,c(P,TIME+1,N))
  store_theta <-array(NA,c(P,TIME+1,N))
  store_alpha <-array(NA,c(PA,TIME+1,N))
  store_aux <-array(NA,c(PA,TIME+1,N))
  #store_omega <- array(NA,c(n,TIME*n,N))
  store_v <-array(NA,c(n+np,TIME,N))
  store_f <-array(NA,c(n,TIME,N))
  store_fc.m<-array(NA,c(n,step_ahead,N))
  store_fc.y<-array(NA,c(n,step_ahead,N))
  store_fc.v<-array(NA,c(n+np,step_ahead,N))
  
  # START GIBBS SAMPLING
  #-------------------------------------------------------------
  
  ptm <- proc.time()
  
  for(it in 1:N){
    
    if(it %% 100==0) cat("\n", it, "loops...") 
    if(it == N) cat("Done!\n") 
    
    # STEP 1: DRAW B
    
    G = phi1*gamma[2:(TIME+1),]
    W = gamma[2:(TIME+1),]*lambda1+(1-gamma[2:(TIME+1),])*lambda0
    
    MATX<-cbind(bigY,G,W,variance,covariance)
    
    mod<-dlm(m0=phi0*gamma[1,],
             C0=diag(gamma[1,]*lambda1/(1-phi1^2)+(1-gamma[1,])*lambda0),
             FF=FF.1,
             GG=GG.1,
             JFF=JFF.1,
             JGG=JGG.1,
             W = W.1,
             JW = JW.1,
             V=V.1,
             JV = JV.1,
             X=MATX)
    
    
    Filtm<-dlmFilter(Y,mod)
    theta<-t(dlmBSample(Filtm))
    
    beta<-matrix(theta,ncol=1,byrow=T)
    betahat<-matrix(theta[,2:(TIME+1)],ncol=1,byrow=T)
    
    thetaind = THETA.B*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA.B*dnorm(beta[1:((TIME+1)*P-P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                      +(1-THETA.B)*dnorm(beta[1:((TIME+1)*P-P),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(beta[1:((TIME+1)*P-P),1]-phi0)
    
    pstar = thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))/(thetaind*dnorm(beta[(P+1):((TIME+1)*P),1],mu,sqrt(lambda1))
                                                                         +(1-thetaind)*dnorm(beta[(P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
    
    pstar[is.na(pstar)]=0
    
    gammino = cpprbinom2(TIME*P,1,pstar)
    
    gammino0 = cpprbinom2(P,1,thetaind[1:P])
    
    gamma = matrix(c(gammino0,gammino),ncol=P,byrow=T)
    
    #gamma[,c(1,1+P/n,1+P/n+P/n)]=1
    
    yhat=oldY-bigX%*%betahat
    
    # STEP 2: DRAW A
    
    hatX = matrix(-yhat,ncol=n,byrow=T)
    hatG = phi1*aux[2:(TIME+1),]
    hatW = aux[2:(TIME+1),]*lambda1+(1-aux[2:(TIME+1),])*lambda0
    Sigma = bigSigma
    
    hat_MATX<-cbind(hatX,hatG,hatW,Sigma)
    
    mod.2<-dlm(m0=phi0*aux[1,],
               C0=diag(aux[1,]*lambda1/(1-phi1^2)+(1-aux[1,])*lambda0,np),
               FF=FF.2,
               GG=GG.2,
               JFF=JFF.2,
               JGG=JGG.2,
               W = W.2,
               JW = JW.2,
               V= V.2,
               JV = JV.2,
               X=hat_MATX)
    
    Filt<-dlmFilter(hatX*10,mod.2)
    a<-t(dlmBSample(Filt))
    
    alpha<-matrix(a,ncol=1,byrow=T)
    alphahat<-matrix(a[,2:(TIME+1)],ncol=1,byrow=T)
    
    auxind = THETA.A*dnorm(alpha[1:((TIME+1)*PA-PA),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA.A*dnorm(alpha[1:((TIME+1)*PA-PA),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                       +(1-THETA.A)*dnorm(alpha[1:((TIME+1)*PA-PA),1],0,sqrt(lambda0)))
    
    mu = phi0+phi1*(alpha[1:((TIME+1)*PA-PA),1]-phi0)
    
    auxstar = auxind*dnorm(alpha[(PA+1):((TIME+1)*PA),1],mu,sqrt(lambda1))/(auxind*dnorm(alpha[(PA+1):((TIME+1)*PA),1],mu,sqrt(lambda1))
                                                                            +(1-auxind)*dnorm(alpha[(PA+1):((TIME+1)*PA),1],0,sqrt(lambda0)))
    auxstar[is.na(auxstar)] = 0
    auxind[is.na(auxind)] = 0
    
    auxgammino = cpprbinom2(TIME*PA,1,auxstar)
    
    auxgammino0 = cpprbinom2(PA,1,auxind[1:PA])
    
    aux = matrix(c(auxgammino0,auxgammino),ncol=PA,byrow=T)
    
    a=a/10
    
    ystar = matrix(NA,ncol=1,nrow=n*TIME)
    FillA<-list()
    giro=1
    for(t in 1:TIME){
      BlockZ = matrix(0, nrow=n, ncol=PA) 
      BlockA = diag(1,nrow=n, ncol=n)
      BlockA[lower.tri(BlockA,diag=F)]=a[,t]
      k=0
      j=1
      l=1
      for(i in 1:(n-1)){
        BlockZ[2+k,j:(j+length(seq(1:i))-1)] = hatX[t,1:i]
        k = k+1
        j = j+length(seq(1:i))
      }
      ystar[giro:(giro+n-1),] = yhat[giro:(giro+n-1),]-BlockZ%*%a[,t]
      giro=giro+n
      FillA[[t]] = solve(BlockA, diag(dim(BlockA)[1]))
    }
    bigA = bdiag_m(FillA)
    
    res = matrix(ystar,ncol=n,byrow=T)
    
    for(i in 1:n){
      res_fast = svsample_fast_cpp(res[,i],
                                   startpara = params, startlatent = start_latent)
      
      
      if(it<20){
        sigtdraw[i,] = rep(-2,dim(res)[1])
        bigSigma[,i] = exp(sigtdraw[i,])
        h.new[,i] = -2
      }else{
        sigtdraw[i,] = res_fast$latent
        bigSigma[,i]=  exp(res_fast$latent)
        sigzeta = res_fast$para[1,3]
        alpha1h = res_fast$para[1,2]
        alpha0h = res_fast$para[1,1]
        h.new[,i] = alpha0h+alpha1h*sigtdraw[i,TIME]+rnorm(1,0,sigzeta)
      }
    }
    
    digSigma = diag(c(matrix(t(bigSigma),nrow=1)))
    
    bigOmega = as.matrix(bigA%*%digSigma%*%t(bigA))
    
    variance = matrix(diag(bigOmega),ncol=n,byrow=T)
    
    if(n==2){
      covariance = bigOmega[row(bigOmega)==col(bigOmega)-1 & col(bigOmega) %% 2==0]
      #cov_firststep = matrix(bigOmega,ncol=n,byrow=T)
      #cov_secondstep = cov_firststep[rowSums(cov_firststep != 0) > 0,]
      #cov_thirdstep = matsplitter(cov_secondstep,n,PA)
      #store_omega[,,it] = t(cov_secondstep)
    }else{
      cov_firststep = matrix(bigOmega,ncol=n,byrow=T)
      cov_secondstep = cov_firststep[rowSums(cov_firststep != 0) > 0,]
      cov_thirdstep = matsplitter(cov_secondstep,n,n)
      covariance = reshape.covariance(cov_thirdstep)
      #store_omega[,,it] = t(cov_secondstep)
    }
    
    
    store_gamma[,,it] = t(gamma) 
    store_theta[,,it] = theta
    store_alpha[,,it] = a
    store_aux[,,it] = t(aux)
    store_v[,,it] = cbind(variance,covariance)
    
    # FORECAST
    
    thetaind.new = THETA.B*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))/(THETA.B*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],phi0,sqrt(lambda1/(1-phi1^2)))
                                                                                                         +(1-THETA.B)*dnorm(beta[(((TIME+1)*P)-P+1):((TIME+1)*P),1],0,sqrt(lambda0)))
    #V.new = bigOmega[(dim(bigOmega)[1]-n+1):dim(bigOmega)[1],(dim(bigOmega)[2]-n+1):dim(bigOmega)[2]]  
    gammino.new = cpprbinom2(P,1,thetaind.new)
    WW.new = diag(gammino.new*lambda1+(1-gammino.new)*lambda0)
    
    mm.new = matrix(Filtm$m[dim(Filtm$m)[1],])
    listC.new <- dlmSvd2var(Filtm$U.C, Filtm$D.C)
    
    variance.new = diag(c(exp(h.new)))
    S.new = diag(c(aux[TIME+1,]*lambda1+(1-aux[TIME+1,])*lambda0),np,np)
    covariance.new = mvrnorm(1,phi1*a[,TIME+1],S.new)
    A.new[lower.tri(A.new)] = covariance.new
    V.new = A.new%*%variance.new%*%t(A.new)
    
    CC.new = listC.new[[length(listC.new)]]
    GG.new = diag(phi1*gammino.new)
    a.new = GG.new%*%mm.new
    f.new = FF.new%*%a.new
    
    R.new = (GG.new)%*%CC.new%*%t(GG.new)+WW.new
    Q.new = FF.new%*%R.new%*%t(FF.new)+V.new
    
    store_fc.m[,,it]<-f.new
    store_fc.v[,,it]<-vechC(Q.new)
    store_fc.y[,,it]<-t(rmvnorm(1,mean=f.new,sigma=Q.new))
    
    
    
  }  
  
  cat("MCMC takes:")
  print(proc.time() - ptm)
  
  return(list(beta=store_theta/10,alpha=store_alpha,
              gamma_beta=store_gamma,gamma_alpha=store_aux,v=store_v,bigSigma=bigSigma,
              fc.m = store_fc.m/10, fc.y = store_fc.y/10, fc.v = store_fc.v/100,
              lags=lags,n=n,N=N,P=P,PA=PA,TIME=TIME))
  
}
