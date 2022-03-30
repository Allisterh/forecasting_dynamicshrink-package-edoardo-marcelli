# Exercise on Airpassengers data
#----------------------------------

# Open library
#-------------------------------
library(dplyr)

# Set working directory 
#-------------------------------

setwd("C:/Users/Edoardo/Dropbox/elementi salvati R")

# Generate Quasi-synthetical data set
#-------------------------------------
sim<-sparse.data.sim(TIME=144,P=4,FP=16,phi0=0,phi1=0.98,lambda1=0.1,v=0,seed=100)
synthetic.y<-sim$y
X<-sim$X
b<-sim$b
y<-log(as.numeric(AirPassengers))*10
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
  assign(paste0('plotairp_',i),Coef.compare.plot(trueb=true_b,column=i,fun=ssvs.bsts,burn=200))
  #paste0('plotairp_',i) %>% save(list=., file = paste0('plotairp_',i,'.Rda'))
}
plot_airpstr = SeasTrendReg.plot(ssvs.bsts,S=11,U=2,burn=200)
save(plot_airpstr,file="plot_airpstr.Rda")
plot(y,type="l")
par(new=T)
plot(seasonality+trend,type="l",col="red")

time=seq(1,144)
df = data.frame(y,trend,seasonality,reg,var.reg,var.s,var.t,time)
plotst501 = ggplot(df,aes(x=time))+
  geom_line(aes(y=y))+
  labs(x="Time",y="")+
  ggtitle("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plotst502 = ggplot(df,aes(x=time))+
  geom_line(aes(y=trend+seasonality+reg),col="red")+
  geom_ribbon(aes(ymin=trend+seasonality+reg-1.96*sqrt(var.s+var.t+var.reg),
                  ymax=trend+seasonality+reg+1.96*sqrt(var.s+var.t+var.reg)),alpha=0.16,fill="red")+
  labs(x="Time",y="")+
  ggtitle("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plotst503 = ggplot(df,aes(x=time))+
  geom_line(aes(y=y))+
  geom_line(aes(y=trend+seasonality+reg),col="red",linetype="dashed")+
  geom_ribbon(aes(ymin=trend+seasonality+reg-1.96*sqrt(var.s+var.t+var.reg),
                  ymax=trend+seasonality+reg+1.96*sqrt(var.s+var.t+var.reg)),alpha=0.16,fill="red")+
  labs(x="Time",y="")+
  ggtitle("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

save(plotst501,file="plotst501.Rda")
save(plotst502,file="plotst502.Rda")
save(plotst503,file="plotst503.Rda")


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

save(plot_compos.y,file="plot_compos.Rda")
save(plot_real.y,file="plot_real.Rda")
save(plot_synthetic.y,file="plot_synthetic.y")