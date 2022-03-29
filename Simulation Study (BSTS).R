# This script produces results of simulations reported in Chapter 3
#-------------------------------------------------------------------

# Open library
#-------------------------------
library(dplyr)

# DO NOT FORGET TO LOAD AUXILIARY FUNCTIONS FROM OTHER SCRIPTS

# Set working directory 
#-------------------------------

setwd("C:/Users/edoar/Dropbox/elementi salvati R")

# Generate Quasi-synthetical data set
#-------------------------------------
sim<-sparse.data.sim(TIME=144,P=4,FP=16,phi0=0,phi1=0.98,lambda1=0.1,v=0,seed=100)
synthetic.y<-sim$y
X<-sim$X
b<-sim$b
real.y<-log(as.numeric(AirPassengers))
y<-real.y*10+synthetic.y #Additive time series decomposition, rescaled by 10
true_b<-sim$true_b
true_ind<-sim$true_ind

# Dynamic EMVS for BSTS model, THETA = 0.2
#------------------------------------------
set.seed(1)
emvs.bsts = DEMVS.Monthly(y,X,N=1000,phi0=0,phi1=0.98,THETA=0.2,
                          lambda0=0.01,lambda1=0.1,
                          sig2u=0.1,sig2d=0.1,sig2tau=0.1,
                          m0d=0,C0d=100,m0u=0,C0u=100,m0tau=0,C0tau=100)
save(emvs.bsts,file="emvsbsts.Rda")
for(i in 1:20){
  assign(paste0('plotbstsem_',i),Coef.compare.plot(trueb = true_b,column = i, fun = emvs.bsts, burn=200, EMVS=T))
  #paste0('plotbstsem_',i) %>% save(list=., file = paste0('plotbstsem_',i,'.Rda'))
}

SeasTrendReg.plot.1(emvs.bsts,frequency="Monthly")

# Dynamic SSVS for BSTS model, THETA = 0.2
#------------------------------------------

set.seed(1)
ssvs.bsts = DSSBSTS_SV(y=y,X=X,N=1000,S=11,U=2,
                       phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                       lambda0=0.01,lambda1=0.1,THETA=0.2)
save(ssvs.bsts,file="ssvsbsts.Rda")
for(i in 1:20){
  assign(paste0('plotbstssv_',i),Coef.compare.plot(trueb=true_b,column=i,fun=ssvs.bsts,burn=200))
  #paste0('plotbstssv_',i) %>% save(list=., file = paste0('plotbstssv_',i,'.Rda'))
}
SeasTrendReg.plot(ssvs.bsts,S=11,U=2,burn=200)

# Full DLM
#-------------------------------------

full.dlm_20 = DSSBSTS_SV(y=y,X=X,N=1000,phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,S=11,U=2,
                         lambda0=0.01,lambda1=0.1,THETA=1)

# Figures
#---------
timeframe = c(1:length(y))
df_qs = data.frame(y,synthetic.y,real.y,timeframe)
plot_synthetic.y <- ggplot(df_qs,aes(x=timeframe))+
  geom_line(aes(y=synthetic.y))+
  labs(x="Time",y="")+
  ggtitle("Synthetic data")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plot_real.y <- ggplot(df_qs,aes(x=timeframe))+
  geom_line(aes(y=real.y*10))+
  labs(x="Time",y="")+
  ggtitle("AirPassengers")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plot_compos.y <- ggplot(df_qs,aes(x=timeframe))+
  geom_line(aes(y=y))+
  labs(x="Time",y="")+
  ggtitle("Resulting Time Series")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

#save(plot_compos.y,file="plot_compos.Rda")
#save(plot_real.y,file="plot_real.Rda")
#save(plot_synthetic.y,file="plot_synthetic.y")

# With 40 Predictors
#----------------------

sim<-sparse.data.sim(TIME=144,P=4,FP=36,phi0=0,phi1=0.98,lambda1=0.1,v=0,seed=100)
synthetic.y<-sim$y
X<-sim$X
b<-sim$b
real.y<-log(as.numeric(AirPassengers))
y<-real.y*10+synthetic.y #Additive time series decomposition, rescaled by 10
true_b_40<-sim$true_b
true_ind_40<-sim$true_ind

emvs.bsts2 = DEMVS.Monthly(y,X,N=1000,phi0=0,phi1=0.98,THETA=0.2,
                           lambda0=0.01,lambda1=0.1,
                           sig2u=0.1,sig2d=0.1,sig2tau=0.1,
                           m0d=0,C0d=100,m0u=0,C0u=100,m0tau=0,C0tau=100)
#save(emvs.bsts2,file="emvsbsts2.Rda")

ssvs.bsts2 = DSSBSTS_SV(y=y,X=X,N=1000,S=11,U=2,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2)
#save(ssvs.bsts2,file="ssvsbsts2.Rda")

# Table
#-------

beta.ssvs_20 = MCMC.out(ssvs.bsts,200)$mean[,1:20]
beta.ssvs_50 = MCMC.out(ssvs.bsts2,200)$mean[,1:40]
noshrink_20 = MCMC.out(full.dlm_20,200)$mean[,1:20]
beta.emvs_20 = t(emvs.bsts$beta[dim(emvs.bsts$beta)[1],,-1])
beta.emvs_50 = t(emvs.bsts2$beta[dim(emvs.bsts2$beta)[1],,-1])
list.bsts = list(noshrink_20,beta.ssvs_20,beta.emvs_20,
                 beta.ssvs_50,beta.emvs_50)

ind_ssvs_20 = MCMC.out(ssvs.bsts,200)$ind[,1:20]
ind_ssvs_50 = MCMC.out(ssvs.bsts2,200)$ind[,1:40]
ind_nosh_20 = MCMC.out(full.dlm_20,200)$ind[,1:20]
ind_emvs_20 = t(emvs.bsts$gamma[dim(emvs.bsts$beta)[1],,-1])
ind_emvs_50 = t(emvs.bsts2$gamma[dim(emvs.bsts2$beta)[1],,-1])
list.ind = list(ind_nosh_20,ind_ssvs_20,ind_emvs_20,
                ind_ssvs_50,ind_emvs_50)

TABLE.3 = matrix(NA,nrow=5,ncol=6)
for(i in c(1:3)){
  TABLE.3[i,1]=SSE(true_b,list.bsts[[i]])
  TABLE.3[i,2]=HmD(true_ind,list.ind[[i]])
  TABLE.3[i,3]=False.Positive(true_ind,list.ind[[i]])
  TABLE.3[i,4]=False.Negative(true_ind,list.ind[[i]])
  TABLE.3[i,5]=False.Active(true_ind,list.ind[[i]])
  TABLE.3[i,6]=False.Non_Active(true_ind,list.ind[[i]])
}
for(i in c(4:5)){
  TABLE.3[i,1]=SSE(true_b_40,list.bsts[[i]])
  TABLE.3[i,2]=HmD(true_ind_40,list.ind[[i]])
  TABLE.3[i,3]=False.Positive(true_ind_40,list.ind[[i]])
  TABLE.3[i,4]=False.Negative(true_ind_40,list.ind[[i]])
  TABLE.3[i,5]=False.Active(true_ind_40,list.ind[[i]])
  TABLE.3[i,6]=False.Non_Active(true_ind_40,list.ind[[i]])
}

TABLE.3 = cbind(matrix(c(20,20,20,40,40)),TABLE.3)
rownames(TABLE.3) = c("Full BSTS (no shrink)","Dynamic SSVS","Dynamic EMVS","Dynamic SSVS", "Dynamic EMVS")
colnames(TABLE.3) = c("p","SSE","Ham. Dist.","FP","FN","FD","FND")
#save(TABLE.3,file="TABLE3.Rda")

# Seasonality, Trend, Slope
#---------------------------

plot_stsl <- SeasTrendSlope.plot(ssvs.bsts,S=11,U=2,burn=200)
#save(plot_stsl,file="plot_stsl.Rda")

