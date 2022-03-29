# This script produces results reported in Chapter 2
#-----------------------------------------------------

# Open library
#-------------------------------
library(dplyr)

# DO NOT FORGET TO LOAD FUNCTIONS FROM THE OTHER SCRIPTS

# Set working directory 
#-------------------------------

setwd("C:/Users/edoar/Dropbox/elementi salvati R")

# Generate a sparse dataset
#--------------------------------

sim<-sparse.data.sim(TIME=144,P=4,FP=46,phi0=0,phi1=0.98,lambda1=0.1,v=0.25,seed=1)

y<-sim$y
X<-sim$X
b<-sim$b
true_b<-sim$true_b
true_ind<-sim$true_ind

# Classic DLM, No Shrinkage, THETA = 1
#--------------------------------------
set.seed(1)
dlm.50<-DSSBSTS_DF(y=y,X=X,N=1000,n0=10,d0=10,
                   phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                   lambda0=0.01,lambda1=0.1,THETA=1,delta=0.9)

for(i in 1:50){
  assign(paste0('plotdlm_',i),Coef.compare.plot(trueb=true_b,column=i,fun=dlm.50,burn=200)) 
  #paste0('plotdlm_',i) %>% save(list=., file = paste0('plotdlm_',i,'.Rda'))
}

# Dynamic SSVS with Disc Fact vol, THETA = 0.2
#----------------------------------------------
set.seed(1)
df.50<-DSSBSTS_DF(y=y,X=X,N=1000,n0=40,d0=10,
                  phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                  lambda0=0.01,lambda1=0.1,THETA=0.2,delta=0.9)

for(i in 1:50){
  assign(paste0('plotdf_',i),Coef.compare.plot(trueb=true_b,column=i,fun=df.50,burn=200)) 
  #paste0('plotdf_',i) %>% save(list=., file = paste0('plotdf_',i,'.Rda'))
}

# Dynamic SSVS with log-AR(1) vol, THETA = 0.2
#----------------------------------------------

set.seed(1)
tvp.reg.sv.50<-DSSBSTS_SV(y=y,X=X,N=1000,
                          phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                          lambda0=0.01,lambda1=0.1,THETA=0.2)

for(i in 1:50){
  assign(paste0('plotsv_',i),Coef.compare.plot(trueb=true_b,column=i,fun=tvp.reg.sv.50,burn=200))
  #paste0('plotsv_',i) %>% save(list=., file = paste0('plotsv_',i,'.Rda'))
}

# Compare volatilites
#---------------------
varplot.df = varplot(df.50,200,lb=0,ub=20)
varplot.sv = varplot(tvp.reg.sv.50,200,lb=0,ub=1)
#save(varplot.sv, file="Varplotsv.Rda")
#save(varplot.df, file="Varplotdf.Rda")

# Dynamic EMVS with disc fact, THETA = 0.2
#------------------------------------------
set.seed(1)
demvs.df = DEMVS(y=y,X=X,N=1000,n0=40,d0=10,
                 phi0=0,phi1=0.98,lambda0=0.01,lambda1=0.1,THETA=0.2,delta=0.9)

for(i in 1:50){
  assign(paste0('plotdemvsdf_',i),Coef.compare.plot(trueb = true_b,column = i, fun = demvs.df, burn=200, EMVS=T))
  #paste0('plotdemvsdf_',i) %>% save(list=., file = paste0('plotdemvsdf_',i,'.Rda'))
}

# Dynamic EMVS with log-N AR(1) vol, THETA = 0.2
#------------------------------------------------
set.seed(1)
demvs.sv = DEMVS_PF(y=y,X=X,N=1000,phi0=0,phi1=0.98,THETA=0.2,lambda0=0.01,lambda1=0.1)

for(i in 1:50){
  assign(paste0('plotdemvssv_',i),Coef.compare.plot(trueb = true_b,column = i, fun = demvs.sv, burn=200, EMVS=T))
  #paste0('plotdemvssv_',i) %>% save(list=., file = paste0('plotdemvssv_',i,'.Rda'))
}

# Volatility comparison
#---------------------
varplot.emvsdf = varplot(demvs.df,200,EMVS=T,lb=0,ub=1)
varplot.emvssv = varplot(demvs.sv,200,EMVS=T,lb=0,ub=1)
#save(varplot.emvssv, file="Varplotemvssv.Rda")
#save(varplot.emvsdf, file="Varplotemvsdf.Rda")

# Table
#-----------------------

betahat_dlm = MCMC.out(dlm.50,200)$mean
betahat_df = MCMC.out(df.50,200)$mean
betahat_sv = MCMC.out(tvp.reg.sv.50,200)$mean
betahat_demdf = t(demvs.df$beta[dim(demvs.df$beta)[1],,-1])
betahat_demsv = t(demvs.sv$beta[dim(demvs.sv$beta)[1],,-1])
list.betahat = list(betahat_dlm,betahat_df,betahat_sv,betahat_demdf,betahat_demsv)

ind_dlm = MCMC.out(dlm.50,200)$ind
ind_df = MCMC.out(df.50,200)$ind
ind_sv = MCMC.out(tvp.reg.sv.50,200)$ind
ind_demdf = t(demvs.df$gamma[dim(demvs.df$beta)[1],,-1])
ind_demsv = t(demvs.sv$gamma[dim(demvs.sv$beta)[1],,-1])
list.ind = list(ind_dlm,ind_df,ind_sv,ind_demdf,ind_demsv)

TABLE.1 = matrix(NA,nrow=5,ncol=6)
for(i in 1:5){
  TABLE.1[i,1]=SSE(trueb,list.betahat[[i]])
  TABLE.1[i,2]=HmD(true_ind,list.ind[[i]])
  TABLE.1[i,3]=False.Positive(true_ind,list.ind[[i]])
  TABLE.1[i,4]=False.Negative(true_ind,list.ind[[i]])
  TABLE.1[i,5]=False.Active(true_ind,list.ind[[i]])
  TABLE.1[i,6]=False.Non_Active(true_ind,list.ind[[i]])
}

colnames(TABLE.1) = c("SSE","Ham. Dist.", "FP","FN","FD","FND") 
rownames(TABLE.1) = c("Full DLM","Dynamic SSVS","Dynamic SSVS (mod.)","Dynamic EMVS","Dynamic EMVS (mod.)")
save(TABLE.1,file="TABLE1.Rda") 

TABLE.2 = matrix(NA,nrow=5,ncol=4)
for(i in 1:5){
  TABLE.2[i,1]=SSE(trueb[,1:4],list.betahat[[i]][,1:4])
  TABLE.2[i,2]=HmD(true_ind[,1:4],list.ind[[i]][,1:4])
  TABLE.2[i,3]=SSE(trueb[,5:50],list.betahat[[i]][,5:50])
  TABLE.2[i,4]=HmD(true_ind[,5:50],list.ind[[i]][,5:50])
}

colnames(TABLE.2) = c("SSE","Ham. Dist.","SSE","Ham. Dist.")
rownames(TABLE.2) = c("Full DLM","Dynamic SSVS","Dynamic SSVS (mod.)","Dynamic EMVS","Dynamic EMVS (mod.)")
#save(TABLE.2,file="TABLE2.Rda")

# Indicator plots

Make_ind_plots = function(y,i,list,true_ind){
  
  i_true = true_ind[,i]
  i_dlm = list.ind[[1]][,i]
  i_df = list.ind[[2]][,i]
  i_sv = list.ind[[3]][,i]
  i_emdf =list.ind[[4]][,i]
  i_emsv = list.ind[[5]][,i]
  
  timeframe = c(1:length(y))
  
  indic_df = data.frame(i_true,i_dlm,i_df,i_sv,i_emdf,i_emsv,timeframe)
  
  ggplot(indic_df,aes(x=timeframe))+
    geom_line(aes(y=i_true), linetype=3,col="black")+
    geom_line(aes(y=i_dlm), linetype=1,col="black")+
    geom_line(aes(y=i_df), linetype=1, col="limegreen")+
    geom_line(aes(y=i_sv), linetype=1, col="darkorange")+
    coord_cartesian(ylim = c(0, +1)) +
    labs(x="Time",
         y="")+
    theme_bw()
}

for(i in 1:6){
  assign(paste0('plotind_',i),Make_ind_plots(y=y,i,list=list.ind,true_ind))
  paste0('plotind_',i) %>% save(list=., file = paste0('plotind_',i,'.Rda'))
}

Make_ind_plots_2 = function(y,i,list,true_ind){
  
  i_true = true_ind[,i]
  i_dlm = list.ind[[1]][,i]
  i_df = list.ind[[2]][,i]
  i_sv = list.ind[[3]][,i]
  i_emdf =list.ind[[4]][,i]
  i_emsv = list.ind[[5]][,i]
  
  timeframe = c(1:length(y))
  
  indic_df = data.frame(i_true,i_dlm,i_df,i_sv,i_emdf,i_emsv,timeframe)
  
  ggplot(indic_df,aes(x=timeframe))+
    geom_line(aes(y=i_true), linetype=3,col="black")+
    geom_line(aes(y=i_emdf), linetype=1, col="limegreen")+
    geom_line(aes(y=i_emsv), linetype=1, col="darkorange")+
    coord_cartesian(ylim = c(0, +1)) +
    labs(x="Time",
         y="")+
    theme_bw()
}

for(i in 1:6){
  assign(paste0('plotind2_',i),Make_ind_plots_2(y=y,i,list=list.ind,true_ind))
  #paste0('plotind2_',i) %>% save(list=., file = paste0('plotind2_',i,'.Rda'))
}

# THE END