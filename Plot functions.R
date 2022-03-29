#-------------------------------------------------------------------#
#                                                                   #
#               Plot functions                                      #
#                                                                   #
#-------------------------------------------------------------------#

# In-sample - Forecast plot
#---------------------------

In_Sample_Forecastplot<-function(y,fun,burn,time){
  
  mcmc.out = MCMC.out(fun,burn)
  
  f = mcmc.out$f
  v = mcmc.out$v
  mean = mcmc.out$mean
  var = mcmc.out$var
  ind = mcmc.out$ind
  q = mcmc.out$q
  varf =mcmc.out$varf
  
  if(missing(time)){
    time = c(1:length(y))
  }else{
    time = time
  }
  
  dataframe = data.frame(y,f,v,q,time,varf)
  
  ggplot(dataframe, aes(x=time))+
    geom_line(aes(y=y))+
    geom_line(aes(y=f),col="red")+
    geom_ribbon(aes(ymin=f-1.96*sqrt(q),ymax=f+1.96*sqrt(q)),alpha=0.16,fill="red")+
    coord_cartesian(ylim = c(-10, +22))+
    labs(x="Time",
         y="")+
    theme_bw()
}

# Var Plot
#---------

varplot<-function(fun,burn,time,lb=0,ub=1,EMVS=T){
  
  if(missing(time)){
    time = c(1:length(y))
  }else{
    time = time
  }
  
  if(missing(EMVS)){
    mcmc.out = MCMC.out(fun,burn)
    
    v = mcmc.out$v
    var.v =mcmc.out$var.v
    
    dataframe = data.frame(v,var.v,time)
    
    ggplot(dataframe, aes(x=time))+
      geom_line(aes(y=sqrt(v)))+
      coord_cartesian(ylim = c(lb, ub))+
      labs(x="Time",
           y="")+
      theme_bw()
    
  }else{
    
    v = fun$v[dim(fun$v)[1],]
    dataframe = data.frame(v,time)
    
    ggplot(dataframe, aes(x=time))+
      geom_line(aes(y=v))+
      coord_cartesian(ylim = c(lb, ub))+
      labs(x="Time",
           y="")+
      theme_bw()
    
  }
}

# Indicator plot
#---------------

Indicatorplot<-function(trueind,fun,burn,time){
  
  mcmc.out = MCMC.out(fun,burn)
  
  f = mcmc.out$f
  v = mcmc.out$v
  mean = mcmc.out$mean
  var = mcmc.out$var
  ind = mcmc.out$ind
  q = mcmc.out$q
  varf =mcmc.out$varf
  
  if(missing(time)){
    time = c(1:length(y))
  }else{
    time = time
  }
  
  dataframe = data.frame(y,trueind,ind,f,v,q,time,varf)
  
  ggplot(dataframe, aes(x=time))+
    geom_line(aes(y=ind[,j]))+
    geom_line(aes(y=trueind[,j]), linetype="dotted")+
    coord_cartesian(ylim = c(0, +1))+
    labs(x="Time",
         y="")+
    theme_bw()
}

# Coefficient plot
#-----------------

Coef.compare.plot<-function(trueb,column,fun,burn,lb,ub,time,EMVS){
  
  require(ggplot2)
  if(missing(time)){
    time = c(1:length(y))
  }else{
    time = time
  }
  
  if(missing(EMVS)){
    mcmc.out = MCMC.out(fun,burn)
    
    P = dim(true_b)[2]
    mean = mcmc.out$mean[,1:P]
    var = mcmc.out$var[,1:P]
    
    j=column
    m = mean[,j]
    C.m = var[,j]
    trueb.m = trueb[,j]
    dataframe = data.frame(m,C.m,time,trueb.m)
    
    ggplot(dataframe, aes(x=time))+
      geom_line(aes(y=trueb.m))+
      geom_line(aes(y=m),col="orange")+
      geom_ribbon(aes(ymin=m-1.96*sqrt(C.m),ymax=m+1.96*sqrt(C.m)),alpha=0.16,fill="orange")+
      coord_cartesian(ylim=c(lb,ub))+
      labs(x="Time",
           y="")+
      theme_bw()
  }else{
    
    j=column
    est.beta = fun$beta
    map_tr = est.beta[dim(est.beta)[1],,]
    m = c(map_tr[j,-1])
    trueb.m = trueb[,j]
    dataframe = data.frame(m,time,trueb.m)
    
    ggplot(dataframe, aes(x=time))+
      geom_line(aes(y=trueb.m))+
      geom_line(aes(y=m),col="green")+
      coord_cartesian(ylim=c(lb,ub))+
      labs(x="Time",
           y="")+
      theme_bw()
    
    
  }
  
}

# Seasonal, Trend, Slope Plot
#-----------------------------

SeasTrendSlope.plot = function(fun,S,U,burn,timeframe){
  
  require(ggpubr)
  require(ggplot2)
  if(missing(timeframe)){
    timeframe = c(1:length(fun$y))
  }else{}
  P=dim(fun$X)[2]
  seasonality = MCMC.out(fun,burn)$mean[,P+1]
  trend = MCMC.out(fun,burn)$mean[,P+S+1]
  slope = MCMC.out(fun,burn)$mean[,P+S+2]
  var.s = MCMC.out(fun,burn)$var[,P+1]
  var.t = MCMC.out(fun,burn)$var[,P+S+1]
  var.sl = MCMC.out(fun,burn)$var[,P+S+2]
  
  sts_df=data.frame(timeframe, seasonality, trend, slope,var.s,var.t,var.sl)
  
  plot_seasonality <-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=seasonality))+
    geom_ribbon(aes(ymin=seasonality-1.96*sqrt(var.s),
                    ymax=seasonality+1.96*sqrt(var.s)),alpha=0.26,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Seasonality")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_trend <-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=trend))+
    geom_ribbon(aes(ymin=trend-1.96*sqrt(var.t),
                    ymax=trend+1.96*sqrt(var.t)),alpha=0.26,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Trend")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_slope<-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=slope))+
    geom_ribbon(aes(ymin=slope-1.96*sqrt(var.sl),
                    ymax=slope+1.96*sqrt(var.sl)),alpha=0.26,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Slope")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggarrange(plot_seasonality,plot_trend,plot_slope,nrow=1)
  
  
  
}

# Seasonal,Trend, Regression Plot (SSVS)
#----------------------------------------

SeasTrendReg.plot = function(fun,S,U,burn,timeframe){
  
  require(ggpubr)
  require(ggplot2)
  
  X = fun$X
  TIME = fun$TIME
  P=dim(fun$X)[2]
  
  if(missing(timeframe)){
    timeframe = c(1:length(fun$y))
  }else{}
  seasonality = MCMC.out(fun,burn)$mean[,P+1]
  trend = MCMC.out(fun,burn)$mean[,P+S+1]
  reg = MCMC.out(fun,burn)$fit.reg
  var.s = MCMC.out(fun,burn)$var[,P+1]
  var.t = MCMC.out(fun,burn)$var[,P+S+1]
  var.reg = MCMC.out(fun,burn)$var.fit.reg
  
  sts_df=data.frame(timeframe, seasonality, trend, reg,var.s,var.t,var.reg)
  
  plot_seasonality <-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=seasonality))+
    geom_ribbon(aes(ymin=seasonality-1.282*sqrt(var.s),
                    ymax=seasonality+1.282*sqrt(var.s)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=seasonality-2.576*sqrt(var.s),
                    ymax=seasonality+2.576*sqrt(var.s)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=seasonality-1.96*sqrt(var.s),
                    ymax=seasonality+1.96*sqrt(var.s)),alpha=0.16,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Seasonality")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_trend <-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=trend))+
    geom_ribbon(aes(ymin=trend-1.96*sqrt(var.t),
                    ymax=trend+1.96*sqrt(var.t)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=trend-1.282*sqrt(var.t),
                    ymax=trend+1.282*sqrt(var.t)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=trend-2.576*sqrt(var.t),
                    ymax=trend+2.576*sqrt(var.t)),alpha=0.16,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Trend")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_reg<-ggplot(sts_df,aes(x=timeframe))+
    geom_line(aes(y=reg))+
    geom_ribbon(aes(ymin=reg-1.96*sqrt(var.reg),
                    ymax=reg+1.96*sqrt(var.reg)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=reg-1.282*sqrt(var.reg),
                    ymax=reg+1.282*sqrt(var.reg)),alpha=0.16,fill="gray40")+
    geom_ribbon(aes(ymin=reg-2.576*sqrt(var.reg),
                    ymax=reg+2.576*sqrt(var.reg)),alpha=0.16,fill="gray40")+
    labs(x="Time",y="")+
    ggtitle("Regression")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggarrange(plot_seasonality,plot_trend,plot_reg,nrow=1)
  
  
  
}

# Seasonal,Trend, Regression Plot (EMVS)
#----------------------------------------

SeasTrendReg.plot.1 = function(fun,frequency,timeframe){
  
  require(ggpubr)
  require(ggplot2)
  
  if(frequency=="Monthly"){
    N=fun$N
    TIME=fun$TIME
    u = fun$u[N,2:dim(fun$u)[2]]
    tau = fun$tau[N,2:dim(fun$u)[2]]
    reg = fun$fit.reg[N,]
    if(missing(timeframe)){
      timeframe = c(1:TIME)
    }else{}
    
    sts_df = data.frame(timeframe,u,tau,reg)
    
    plot_seasonality <-ggplot(sts_df,aes(x=timeframe))+
      geom_line(aes(y=tau))+
      labs(x="Time",y="")+
      ggtitle("Seasonality")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    plot_trend <-ggplot(sts_df,aes(x=timeframe))+
      geom_line(aes(y=u))+
      labs(x="Time",y="")+
      ggtitle("Trend")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    plot_reg<-ggplot(sts_df,aes(x=timeframe))+
      geom_line(aes(y=reg))+
      labs(x="Time",y="")+
      ggtitle("Regression")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ggarrange(plot_seasonality,plot_trend,plot_reg,nrow=1)
  }else{}
  
  
  
}

# How many variables activate? Check this with this function
#------------------------------------------------------------

Active.coef.plot = function(fun, burn, EMVS=F,timeframe){
  
  if(missing(timeframe)){
    TIME=fun$TIME
    timeframe = c(1:TIME)
  }else{}
  
  if(EMVS==F){
    mcmc = MCMC.out(fun, burn)
    aux_matrix = mcmc$ind
    aux_matrix[aux_matrix>0.5]=1
    aux_matrix[aux_matrix<0.5]=0
    aux_vec = rowSums(aux_matrix)
    aux_df = data.frame(aux_vec,timeframe)
  }else{
    
    aux_matrix=fun$gamma[fun$N,,] 
    aux_matrix[aux_matrix>0.5]=1
    aux_matrix[aux_matrix<0.5]=0
    aux_vec = rowSums(t(aux_matrix)[-1,])
    aux_df = data.frame(aux_vec,timeframe)
    
    
  }
  
  ggplot(aux_df,aes(x=timeframe))+
    geom_line(aes(y=aux_vec))+
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    labs(x="Time",y="")+
    coord_cartesian(ylim = c(0, max(aux_vec)+1))+
    theme_bw()
}

# Forecast function (for univariate time series models)
#------------------------------------------------------

Forecastplot = function(y,from,to,start_fc,fc.f,fc.Q,ub,lb,timeframe){
  
  train = start_fc-from
  test = to-start_fc
  
  if(missing(timeframe)){
    timeframe = c(1:length(y))
  }else{}
  
  if(missing(fc.Q)){
    
    vec.fcst=c(rep(NA,train),fc.f)
    vec.y_before = c(y[from:start_fc],rep(NA,test))
    vec.y_after = c(rep(NA,train),y[start_fc:to])
    lil.timeframe=timeframe[from:to]
    
    fcst_df = data.frame(vec.y_before,vec.y_after,vec.fcst,lil.timeframe)
    
    
    ggplot(fcst_df,aes(x=lil.timeframe))+
      geom_line(aes(y=vec.y_before))+
      geom_line(aes(y=vec.y_after),linetype="dashed")+
      geom_line(aes(y=vec.fcst),col="blue")+labs(x="Time",
                                                 y="")+
      theme_bw()
    
  }else{
    
    vec.fcst=c(rep(NA,train),fc.f)
    qvec.fcst =c(rep(NA,train),fc.Q)
    vec.y_before = c(y[from:start_fc],rep(NA,test))
    vec.y_after = c(rep(NA,train),y[start_fc:to])
    lil.timeframe=timeframe[from:to]
    
    fcst_df = data.frame(vec.y_before,vec.y_after,vec.fcst,qvec.fcst,lil.timeframe)
    
    ggplot(fcst_df,aes(x=lil.timeframe))+
      geom_line(aes(y=vec.y_before))+
      geom_line(aes(y=vec.y_after),linetype="dashed")+
      geom_line(aes(y=vec.fcst),col="blue")+
      geom_ribbon(aes(ymin=vec.fcst-1.96*sqrt(qvec.fcst),
                      ymax=vec.fcst+1.96*sqrt(qvec.fcst)),alpha=0.16,fill="blue")+
      geom_ribbon(aes(ymin=vec.fcst-1.282*sqrt(qvec.fcst),
                      ymax=vec.fcst+1.282*sqrt(qvec.fcst)),alpha=0.16,fill="blue")+
      geom_ribbon(aes(ymin=vec.fcst-2.576*sqrt(qvec.fcst),
                      ymax=vec.fcst+2.576*sqrt(qvec.fcst)),alpha=0.16,fill="blue")+
      coord_cartesian(ylim = c(lb, ub))+
      labs(x="Time",
           y="")+
      theme_bw()
  }
  
}

# Forecast function (for multivariate time series models)
#--------------------------------------------------------

Forecastplot.TVP.VAR = function(from,to,start_fc,X,fc_m,fc_v,i,lb,ub,title,timeframe){
  require(ggplot2)
  train = start_fc-from
  test = to-start_fc
  
  if(missing(timeframe)){
    timeframe=c(1:length(from:to))}else{
      timeframe=timeframe[from:to]
    }
  
  variance_t = c()
  for(t in 1:length(start_fc:to)){
    v = matrix(fc_v[,t])
    mat=matrix(NA,3,3)
    mat[lower.tri(mat,diag=T)]=v
    mat[upper.tri(mat)]=mat[lower.tri(mat)]
    diag_mat = diag(mat)
    variance_t[t] = diag_mat[i]
  }
  
  variance = c(rep(NA,train),variance_t)
  mean = c(rep(NA,train),fc_m[i,])
  obs_f = c(rep(NA,train),X[start_fc:to,i])
  obs_p = c(X[from:start_fc,i],rep(NA,test))
  
  fcst_df = data.frame(timeframe,variance,mean,obs_f,obs_p)
  
  ggplot(fcst_df,aes(x=timeframe))+
    geom_line(aes(y=obs_p))+
    geom_line(aes(y=obs_f),linetype="dashed")+
    geom_line(aes(y=mean),col="blue")+
    geom_ribbon(aes(ymin=mean-1.96*sqrt(variance),
                    ymax=mean+1.96*sqrt(variance)),alpha=0.16,fill="blue")+
    geom_ribbon(aes(ymin=mean-1.282*sqrt(variance),
                    ymax=mean+1.282*sqrt(variance)),alpha=0.16,fill="blue")+
    geom_ribbon(aes(ymin=mean-2.576*sqrt(variance),
                    ymax=mean+2.576*sqrt(variance)),alpha=0.16,fill="blue")+
    ggtitle(title)+
    coord_cartesian(ylim = c(lb, ub))+
    labs(x="Time",
         y="")+
    theme_bw()
  
}