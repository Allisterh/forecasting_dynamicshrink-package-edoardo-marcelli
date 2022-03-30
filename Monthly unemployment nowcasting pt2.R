# Monthly Unemployment Forecasting
#----------------------------------
library(tidyverse)
library(fbi)
library(fredr)
library(BatchGetSymbols)
library(quantmod)
setwd("C:/Users/Edoardo/Dropbox/elementi salvati R")

# Data
#----------------------------------------------------------------------------

Monthly_GoogleTrends<-read.csv("Monthly_GoogleTrends.csv")
monthly_un_nsa<-read.csv("UNRATENSA_Monthly_NSA.csv")
fredr_set_key("583ff02e600f88e535efeec886cb801a")

nikkei_fred = fredr(
  series_id = "NIKKEI225",
  observation_start = as.Date("2004-01-01"),
  observation_end = as.Date("2021-01-01"),
  frequency = "m", # quarterly
  units = "lin" # change over previous value
)
nikkei = nikkei_fred$value

nasdaq_fred = fredr(
  series_id = "NASDAQ100",
  observation_start = as.Date("2004-01-01"),
  observation_end = as.Date("2021-01-01"),
  frequency = "m", # quarterly
  units = "lin" # change over previous value
)
nasdaq = nasdaq_fred$value

DJI_yahoo = getSymbols(Symbols = '^DJI', auto.assign = FALSE, src = "yahoo",
                       from = "2004-01-01", to = "2021-02-01",periodicity = "monthly")

DJI = as.vector(DJI_yahoo$DJI.Close)

#-----------------------------------------------------------------------------
# Dependent variable

LAG0<-monthly_un_nsa[monthly_un_nsa$DATE >= "2004-01-01" & monthly_un_nsa$DATE <= "2021-01-01", ]
y<-as.numeric(LAG0$UNRATENSA)
y = scale(y)*10
timeframe = as.Date(LAG0$DATE)

# Lagged values (y_t-1)

LAG1 = monthly_un_nsa[monthly_un_nsa$DATE >= "2003-12-01" & monthly_un_nsa$DATE <= "2020-12-01", ]
y_t_1 = as.numeric(LAG1$UNRATENSA)

LAG2 = monthly_un_nsa[monthly_un_nsa$DATE >= "2003-09-01" & monthly_un_nsa$DATE <= "2020-09-01", ]
y_t_2 = as.numeric(LAG1$UNRATENSA)

# Google Trends Regressors (Z_t)

Z_t = Monthly_GoogleTrends[Monthly_GoogleTrends$out.date >= "2004-01-01" & Monthly_GoogleTrends$out.date <= "2021-01-01", ]
Z_t = Z_t %>% dplyr::select(-c(X,out.date))

# Potential Leading Indicators

data_fredmd = fredmd(file = "current.csv", date_start = NULL, date_end = NULL, transform = F)

# Contemporaneus Regressors (W_t)

contemporaneus_reg = data_fredmd[data_fredmd$date >= "2004-01-01" & data_fredmd$date <= "2021-01-01", ]
Cont.reg = data.frame(
  contemporaneus_reg$`S&P 500`,
  contemporaneus_reg$OILPRICEx
  #contemporaneus_reg$CLAIMSx
)

Cont.reg.1 = as.SASMAT(Cont.reg)
W_t = cbind(Cont.reg.1,nikkei,DJI,nasdaq)

# Potential Leading Indicators

lagged_reg = data_fredmd[data_fredmd$date >= "2003-12-01" & data_fredmd$date <= "2020-12-01", ]
Lead.Ind = data.frame(
  lagged_reg$CPIAUCSL,
  lagged_reg$INDPRO,
  lagged_reg$TWEXAFEGSMTHx,
  lagged_reg$RETAILx,
  lagged_reg$DPCERA3M086SBEA,
  lagged_reg$CMRMTSPLx,
  lagged_reg$HOUST,
  lagged_reg$PAYEMS,
  lagged_reg$RPI
)

X_t_1 = as.matrix(Lead.Ind)
# Final Matrix X
sub_Z_t = Z_t[,c("unemployment",
                 "federal.unemployment",
                 "compensation.package",
                 "Unemployment.compensation",
                 "Unemployment.agency",
                 "employee.benefits",
                 "unemployment.check",
                 "unemployment.statistics",
                 "unemployment.pa",
                 "unemployment.office",
                 "unemployment.insurance",
                 "unemployment.depression",
                 "unemployment.benefits",
                 "subsidies")]

#-------------------------------------------------------------------------
# End dataset building
#-------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Seasonal Adjustment (some series do not adjust, please skip those ones)
#--------------------------------------------------------------------------


sub_X = cbind(y_t_1,W_t,X_t_1,sub_Z_t)
matrix = sub_X
require(seasonal)
matrix_sa=matrix(NA,nrow=dim(matrix)[1],ncol=dim(matrix)[2])
for(i in 1:dim(matrix)[2]){
  print(i)
  seas_adj <- seas(ts(matrix[,i],start=2000,frequency=12),x11 = "")
  matrix_sa[,i] = as.numeric(final(seas_adj))
}
X_sa = matrix_sa[, colSums(is.na(matrix_sa)) != nrow(matrix_sa)]
X = X_sa
X = data.matrix(X)
X = scale(X)

# We are ready to analysis

#---------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------

# Dynamic EMVS full sample

unemp_demvs = DEMVS.Monthly(y=y,X=X,N=1000,phi0=0,phi1=0.98,lambda0=0.01,lambda1=0.1,THETA=0.2)
#save(unemp_demvs,file="unemp_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/unemp_demvs.Rda")

# how many coefficient activates?

for(i in 1:dim(X)[2]){
  if(max(unemp_demvs$gamma[1000,i,])>0.50){print(i)}else{}
}
STRplot_demvs = SeasTrendReg.plot.1(unemp_demvs,"Monthly",timeframe)
#save(STRplot_demvs,file="STRplot_demvs.Rda")
ACP_demvs = Active.coef.plot(unemp_demvs,EMVS=T,timeframe)
#save(ACP_demvs,file="ACP_demvs.Rda")

# Dynamic SSVS full sample

unemp_dssbsts = DSSBSTS_SV(y=y,X=X,N=1000,U=2,S=11,
                           phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                           lambda0=0.01,lambda1=0.1,THETA=0.2)

#save(unemp_dssbsts,file="unemp_dssbsts.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/unemp_dssbsts.Rda")

VP_bsts = varplot(unemp_dssbsts,burn=200,time=timeframe,lb=0,ub=0.1)
#save(VP_bsts,file="VP_bsts.Rda")

STRplot_bsts = SeasTrendReg.plot(unemp_dssbsts,U=2,S=11,burn=200)
#save(STRplot_bsts,file="STRplot_bsts.Rda")
ACP_bsts = Active.coef.plot(unemp_dssbsts,burn=200,timeframe=timeframe)
#save(ACP_bsts,file="ACP_bsts.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/STRplot_bsts.Rda")

# how many coefficient activates?

mcmc = MCMC.out(unemp_dssbsts,200)

for(i in 1:(dim(X)[2])){
  if(max(mcmc$ind[,i])>0.50){print(i)}else{}
}

#-----------------------------------------------------------------------
# Forecasting 
#-----------------------------------------------------------------------

# from first half of the sample

# THETA = 0.2

sas.theta02<-c()
qsas.theta02<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=200,U=2,S=11,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas.theta02[j]<-mean(first.try$fc.f)
  qsas.theta02[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
save(sas.theta02,file="sas_theta02.Rda")
save(qsas.theta02,file="qsas_theta02.Rda")
load("C:/Users/Edoardo/Dropbox/elementi salvati R/sas_theta02.Rda")
load("C:/Users/Edoardo/Dropbox/elementi salvati R/qsas_theta02.Rda")
mase1 = mase(y[96:116],sas.theta02)
wmape1 = wmape(y[96:116],sas.theta02)
rmse1 = rmse(y[96:116],sas.theta02)
mae1 = mae(y[96:116],sas.theta02)


fcst_plot_1 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,sas.theta02,qsas.theta02,timeframe)
#save(fcst_plot_1,file="fcst_plot_1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/fcst_plot_1.Rda")



lpds(y[95:115],sas.theta02_ns,qsas.theta02_ns)
lpds(y[96:116],sas.theta02,qsas.theta02)

# THETA 1

sas.theta1<-c()
qsas.theta1<-c()
j=1
for(i in 95:115){
  print(i)
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas.theta1[j]<-mean(second.try$fc.f)
  qsas.theta1[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(sas.theta1,file="sas_theta1.Rda")
#save(qsas.theta1,file="qsas_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qsas_theta1.Rda")
mase2 = mase(y[96:116],sas.theta1)
wmape2 = wmape(y[96:116],sas.theta1)
rmse2 = rmse(y[96:116],sas.theta1)
mae2 = mae(y[96:116],sas.theta1)

fcst_plot_2 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,sas.theta1,qsas.theta1,timeframe)
#save(fcst_plot_2,file="fcst_plot_2.Rda")

# THETA 0.2

sas.theta02_demvs<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           sig2u=1,sig2d=1,sig2tau=1,
                           lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas.theta02_demvs[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(sas.theta02_demvs,file="sas_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta02_demvs.Rda")
mase3 = mase(y[96:116],sas.theta02_demvs)
wmape3 = wmape(y[96:116],sas.theta02_demvs)
rmse3 = rmse(y[96:116],sas.theta02_demvs)
mae3 = mae(y[96:116],sas.theta02_demvs)

fcst_plot_3 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,fc.f = sas.theta02_demvs,timeframe = timeframe)
#save(fcst_plot_3,file="fcst_plot_3.Rda")


# THETA 1

sas.theta1_demvs<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas.theta1_demvs[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
save(sas.theta1_demvs,file="sas_theta1_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1_demvs.Rda")
mase3 = mase(y[96:116],sas.theta1_demvs)
wmape3 = wmape(y[96:116],sas.theta1_demvs)
rmse3 = rmse(y[96:116],sas.theta1_demvs)
mae3 = mae(y[96:116],sas.theta1_demvs)

fcst_plot_4 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,fc.f = sas.theta1_demvs,timeframe = timeframe)
#save(fcst_plot_4,file="fcst_plot_4.Rda")

# TVP-reg
#--------------------------------------------------------

X = cbind(matrix(1,ncol=1,nrow=dim(X)[1]),X)

# THETA = 0.2

sas.theta02_ns<-c()
qsas.theta02_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X,activate=T)
  sas.theta02_ns[j]<-mean(first.try$fc.f)
  qsas.theta02_ns[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
#save(sas.theta02_ns,file="sas_theta02_ns.Rda")
#save(qsas.theta02_ns,file="qsas_theta02_ns.Rda")
mase4 = mase(y[96:116],sas.theta02_ns)
wmape4 = wmape(y[96:116],sas.theta02_ns)
rmse4 = rmse(y[96:116],sas.theta02_ns)
mae4 = mae(y[96:116],sas.theta02_ns)
load("qsas_theta1_ns.Rda")


fcst_plot_5 = Forecastplot(y,from=86,to=116,start_fc=96,ub=25,lb=0,
                           fc.f = sas.theta02_ns,fc.Q = qsas.theta02_ns,timeframe = timeframe)
#save(fcst_plot_5,file="fcst_plot_5.Rda")

# THETA 1

sas.theta1_ns<-c()
qsas.theta1_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas.theta1_ns[j]<-mean(second.try$fc.f)
  qsas.theta1_ns[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(sas.theta1_ns,file="sas_theta1_ns.Rda")
#save(qsas.theta1_ns,file="qsas_theta1_ns.Rda")
mase5 = mase(y[96:116],sas.theta1_ns)
wmape5 = wmape(y[96:116],sas.theta1_ns)
rmse5 = rmse(y[96:116],sas.theta1_ns)
mae5 = mae(y[96:116],sas.theta1_ns)

fcst_plot_6 = Forecastplot(y,from=86,to=116,start_fc=96,ub=25,lb=0,
                           fc.f = sas.theta1_ns,fc.Q = qsas.theta1_ns,
                           timeframe = timeframe)
#save(fcst_plot_6,file="fcst_plot_6.Rda")

# THETA 0.2

sas.theta02_demvs_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas.theta02_demvs_ns[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(sas.theta02_demvs_ns,file="sas_theta02_demvs_ns.Rda")
mase6 = mase(y[96:116],sas.theta02_demvs_ns)
wmape6 = wmape(y[96:116],sas.theta02_demvs_ns)
rmse6 = rmse(y[96:116],sas.theta02_demvs_ns)
mae6 = mae(y[96:116],sas.theta02_demvs_ns)

fcst_plot_7 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,
                           fc.f=sas.theta02_demvs_ns,timeframe=timeframe)
#save(fcst_plot_7,file="fcst_plot_7.Rda")


# THETA 1

sas.theta1_demvs_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas.theta1_demvs_ns[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
#save(sas.theta1_demvs_ns,file="sas_theta1_demvs_ns.Rda")
mase7 = mase(y[96:116],sas.theta1_demvs_ns)
wmape7 = wmape(y[96:116],sas.theta1_demvs_ns)
rmse7 = rmse(y[96:116],sas.theta1_demvs_ns)
mae7 = mae(y[96:116],sas.theta1_demvs_ns)

fcst_plot_8 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,
                           fc.f=sas.theta1_demvs_ns,timeframe=timeframe)
#save(fcst_plot_8,file="fcst_plot_8.Rda")


#...................................................................
# Break
#...................................................................

# THETA = 0.2

sas2.theta02<-c()
qsas2.theta02<-c()
j=1
for(i in 183:203){
  
  old_y<-y[100:i]
  old_X<-X[100:i,]
  new_X<-X[i+1,]
  #first.try = BSTS_SeasTrend(y,S=11,U=2,sig2u=0.1,sig2d=0.1,sig2tau=0.1,N=100,n0=10,d0=10,delta=0.9,v0=0.25)
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas2.theta02[j]<-mean(first.try$fc.f)
  qsas2.theta02[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
save(sas2.theta02,file="sas2_theta02.Rda")
save(qsas2.theta02,file="qsas2_theta02.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta02.Rda")
mase11 = mase(y[184:204],sas2.theta02)
wmape11 = wmape(y[184:204],sas2.theta02)
rmse11 = rmse(y[184:204],sas2.theta02)
mae11 = mae(y[184:204],sas2.theta02)


fcst_plot_11 = Forecastplot(y,from=174,to=204,start_fc=184,lb=-15,ub=50,sas2.theta02,qsas2.theta02,timeframe)
#save(fcst_plot_11,file="fcst_plot_11.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/fcst_plot_11.Rda")

# THETA 1

sas2.theta1<-c()
qsas2.theta1<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas2.theta1[j]<-mean(second.try$fc.f)
  qsas2.theta1[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
save(sas2.theta1,file="sas2_theta1.Rda")
save(qsas2.theta1,file="qsas2_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta1.Rda")
mase12 = mase(y[184:204],sas2.theta1)
wmape12 = wmape(y[184:204],sas2.theta1)
rmse12 = rmse(y[184:204],sas2.theta1)
mae12 = mae(y[184:204],sas2.theta1)

fcst_plot_12 = Forecastplot(y,from=174,to=204,start_fc=184,lb=-15,ub=50,sas2.theta1,qsas2.theta1,timeframe)
#save(fcst_plot_12,file="fcst_plot_12.Rda")

# THETA 0.2

sas2.theta02_demvs<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas2.theta02_demvs[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
save(sas2.theta02_demvs,file="sas2_theta02_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02_demvs.Rda")
mase13 = mase(y[184:204],sas2.theta02_demvs)
wmape13 = wmape(y[184:204],sas2.theta02_demvs)
rmse13 = rmse(y[184:204],sas2.theta02_demvs)
mae13 = mae(y[184:204],sas2.theta02_demvs)

fcst_plot_13 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,fc.f = sas2.theta02_demvs,timeframe = timeframe)
#save(fcst_plot_13,file="fcst_plot_13.Rda")


# THETA 1

sas2.theta1_demvs<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas2.theta1_demvs[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
save(sas2.theta1_demvs,file="sas2_theta1_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1_demvs.Rda")
mase14 = mase(y[184:204],sas2.theta1_demvs)
wmape14 = wmape(y[184:204],sas2.theta1_demvs)
rmse14 = rmse(y[184:204],sas2.theta1_demvs)
mae14 = mae(y[184:204],sas2.theta1_demvs)

fcst_plot_14 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,fc.f = sas2.theta1_demvs,timeframe = timeframe)
#save(fcst_plot_14,file="fcst_plot_14.Rda")

# TVP-reg
#--------------------------------------------------------

X = cbind(matrix(1,ncol=1,nrow=dim(X)[1]),X)

# THETA = 0.2

sas2.theta02_ns<-c()
qsas2.theta02_ns<-c()
j=1
for(i in  183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X,activate=F)
  sas2.theta02_ns[j]<-mean(first.try$fc.f)
  qsas2.theta02_ns[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
save(sas2.theta02_ns,file="sas2_theta02_ns.Rda")
save(qsas2.theta02_ns,file="qsas2_theta02_ns.Rda")
mase15 = mase(y[184:204],sas2.theta02_ns)
wmape15 = wmape(y[184:204],sas2.theta02_ns)
rmse15 = rmse(y[184:204],sas2.theta02_ns)
mae15 = mae(y[184:204],sas2.theta02_ns)


fcst_plot_15 = Forecastplot(y,from=174,to=204,start_fc=184,ub=50,lb=-15,
                            fc.f = sas2.theta02_ns,fc.Q = qsas2.theta02_ns,timeframe = timeframe)
#save(fcst_plot_15,file="fcst_plot_15.Rda")

# THETA 1

sas2.theta1_ns<-c()
qsas2.theta1_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas2.theta1_ns[j]<-mean(second.try$fc.f)
  qsas2.theta1_ns[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
save(sas2.theta1_ns,file="sas2_theta1_ns.Rda")
save(qsas2.theta1_ns,file="qsas2_theta1_ns.Rda")
mase16 = mase(y[184:204],sas2.theta1_ns)
wmape16 = wmape(y[184:204],sas2.theta1_ns)
rmse16 = rmse(y[184:204],sas2.theta1_ns)
mae16 = mae(y[184:204],sas2.theta1_ns)

fcst_plot_16 = Forecastplot(y,from=174,to=204,start_fc=184,ub=50,lb=-10,
                            fc.f = sas2.theta1_ns,fc.Q = qsas2.theta1_ns,
                            timeframe = timeframe)
#save(fcst_plot_16,file="fcst_plot_16.Rda")

# THETA 0.2

sas2.theta02_demvs_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  sas2.theta02_demvs_ns[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(sas2.theta02_demvs_ns,file="sas2_theta02_demvs_ns.Rda")
mase17 = mase(y[184:204],sas2.theta02_demvs_ns)
wmape17 = wmape(y[184:204],sas2.theta02_demvs_ns)
rmse17 = rmse(y[184:204],sas2.theta02_demvs_ns)
mae17 = mae(y[184:204],sas2.theta02_demvs_ns)

fcst_plot_17 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,
                            fc.f=sas2.theta02_demvs_ns,timeframe=timeframe)
#save(fcst_plot_17,file="fcst_plot_7.Rda")


# THETA 1

sas2.theta1_demvs_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  sas2.theta1_demvs_ns[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}

save(sas2.theta1_demvs_ns,file="sas2_theta1_demvs_ns.Rda")
mase18 = mase(y[184:204],sas2.theta1_demvs_ns)
wmape18 = wmape(y[184:204],sas2.theta1_demvs_ns)
rmse18 = rmse(y[184:204],sas2.theta1_demvs_ns)
mae18 = mae(y[184:204],sas2.theta1_demvs_ns)

fcst_plot_18 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,
                            fc.f=sas2.theta1_demvs_ns,timeframe=timeframe)
#save(fcst_plot_18,file="fcst_plot_8.Rda")
Active.coef.plot(forth.try,200,EMVS=T,timeframe=timeframe)

#--------------------------------------------------------------------
# End Forecasting
#--------------------------------------------------------------------

# Break

#--------------------------------------------------------------------
# Generate Table
#-------------------------------------------------------------------

load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qsas2_theta1_ns.Rda")

load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta02_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas_theta1_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta02_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/sas2_theta1_demvs_ns.Rda")



list_sas = list(sas.theta02,sas.theta1,sas.theta02_ns,sas.theta1_ns,
                sas.theta02_demvs,sas.theta1_demvs,sas.theta02_demvs_ns,
                sas.theta1_demvs_ns)

list_qsas = list(qsas.theta02,qsas.theta1,qsas.theta02_ns,qsas.theta1_ns)

list_sas2 = list(sas2.theta02,sas2.theta1,sas2.theta02_ns,sas2.theta1_ns,
                 sas2.theta02_demvs,sas2.theta1_demvs,sas2.theta02_demvs_ns,
                 sas2.theta1_demvs_ns)

list_qsas2 = list(qsas2.theta02,qsas2.theta1,qsas2.theta02_ns,qsas2.theta1_ns)

SASMAT = matrix(NA,ncol = 8,nrow = 8)
n = 1
for(i in list_sas){
  
  SASMAT[n,1] = rmse(y[96:116],i)
  SASMAT[n,2] = wmape(y[96:116],i)
  SASMAT[n,3] = mase(y[96:116],i)
  n=n+1
  
}

SASMAT[1,4] = lpds(y[96:116],sas.theta02,qsas.theta02)
SASMAT[2,4] = lpds(y[96:116],sas.theta1,qsas.theta1)
SASMAT[3,4] = lpds(y[96:116],sas.theta02_ns,qsas.theta02_ns)
SASMAT[4,4] = lpds(y[96:116],sas.theta1_ns,qsas.theta1_ns)


n = 1
for(i in list_sas2){
  
  SASMAT[n,5] = rmse(y[184:204],i)
  SASMAT[n,6] = wmape(y[184:204],i)
  SASMAT[n,7] = mase(y[184:204],i)
  n=n+1
  
}

SASMAT[1,8] = lpds(y[184:204],sas2.theta02,qsas2.theta02)
SASMAT[2,8] = lpds(y[184:204],sas2.theta1,qsas2.theta1)
SASMAT[3,8] = lpds(y[184:204],sas2.theta02_ns,qsas2.theta02_ns)
SASMAT[4,8] = lpds(y[184:204],sas2.theta1_ns,qsas2.theta1_ns)

save(SASMAT,file="SASMAT.Rda")