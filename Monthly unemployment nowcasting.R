# Monthly Unemployment Forecasting
#----------------------------------------------------------------------------

library(tidyverse)
library(fbi)
library(fredr)
library(BatchGetSymbols)
library(quantmod)
setwd("C:/Users/edoar/Dropbox/elementi salvati R")

#----------------------------------------------------------------------------
# Dataset building
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
)

Cont.reg.1 = as.matrix(Cont.reg)
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

list.predictors = c("Dow Jones","NIKKEI225","NASDAQ100",
                    "SP500","OILPRICEx",
                    "CPIAUSL","INDPRO","TWEXAFEGSMTHx","RETAILx",
                    "DPCERA3M086SBEA","CMRMTSPLx","HOUST","PAYEMS","RPI",
                    "unemployment",
                    "federal unemployment",
                    "compensation package",
                    "Unemployment compensation",
                    "Unemployment agency",
                    "employee benefits",
                    "unemployment check",
                    "unemployment statistics",
                    "unemployment pa",
                    "unemployment office",
                    "unemployment insurance",
                    "unemployment depression",
                    "unemployment benefits",
                    "subsidies")

List.Predictors = matrix(list.predictors)
#save(List.Predictors,file="List_predictors.Rda")

sub_X = cbind(y_t_1,W_t,X_t_1,sub_Z_t)
X = data.matrix(sub_X)
X = scale(X)

#---------------------------------------------------------------------------
# End dataset building
#---------------------------------------------------------------------------

# We are ready for the analysis

#---------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------

# Dynamic EMVS full sample

unemp_demvs = DEMVS.Monthly(y=y,X=X,N=1000,phi0=0,phi1=0.98,
                            lambda0=0.01,lambda1=0.1,THETA=0.2)
#save(unemp_demvs,file="unemp_demvs.Rda")

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

VP_bsts = varplot(unemp_dssbsts,burn=200,time=timeframe,lb=0,ub=1)
#save(VP_bsts,file="VP_bsts.Rda")

STRplot_bsts = SeasTrendReg.plot(unemp_dssbsts,U=2,S=11,burn=200,timeframe)
#save(STRplot_bsts,file="STRplot_bsts.Rda")
ACP_bsts = Active.coef.plot(unemp_dssbsts,burn=200,timeframe=timeframe)
#save(ACP_bsts,file="ACP_bsts.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/STRplot_bsts.Rda")

# how many coefficient activates? #Answ:16,17,22,24,25,27

mcmc = MCMC.out(unemp_dssbsts,200)

for(i in 1:(dim(X)[2])){
  if(max(mcmc$ind[,i])>0.50){print(i)}else{}
}

#-----------------------------------------------------------------------
# End of the analysis
#-----------------------------------------------------------------------

# Let us step into forecasting

#----------------------------------------------------------------------
# Forecasting 
#----------------------------------------------------------------------

# from first half of the sample

# THETA = 0.2

vec.theta02<-c()
qvec.theta02<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec.theta02[j]<-mean(first.try$fc.f)
  qvec.theta02[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
#save(vec.theta02,file="vec_theta02.Rda")
#save(qvec.theta02,file="qvec_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta02.Rda")
mase(y[96:116],vec.theta02)
wmape(y[96:116],vec.theta02)
rmse(y[96:116],vec.theta02)
mae(y[96:116],vec.theta02)


fcst_plot_1 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,vec.theta02,qvec.theta02,timeframe)
#save(fcst_plot_1,file="fcst_plot_1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/fcst_plot_1.Rda")

lpds(y[95:115],vec.theta02_ns,qvec.theta02_ns)
lpds(y[96:116],vec.theta02,qvec.theta02)

# THETA 1

vec.theta1<-c()
qvec.theta1<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec.theta1[j]<-mean(second.try$fc.f)
  qvec.theta1[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(vec.theta1,file="vec_theta1.Rda")
#save(qvec.theta1,file="qvec_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta1.Rda")
mase(y[96:116],vec.theta1)
wmape(y[96:116],vec.theta1)
rmse(y[96:116],vec.theta1)
mae(y[96:116],vec.theta1)

fcst_plot_2 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,vec.theta1,qvec.theta1,timeframe)
#save(fcst_plot_2,file="fcst_plot_2.Rda")

# THETA 0.2

vec.theta02_demvs<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec.theta02_demvs[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(vec.theta02_demvs,file="vec_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02_demvs.Rda")
mase3 = mase(y[96:116],vec.theta02_demvs)
wmape3 = wmape(y[96:116],vec.theta02_demvs)
rmse3 = rmse(y[96:116],vec.theta02_demvs)
mae3 = mae(y[96:116],vec.theta02_demvs)

fcst_plot_3 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,fc.f = vec.theta02_demvs,timeframe = timeframe)
#save(fcst_plot_3,file="fcst_plot_3.Rda")


# THETA 1

vec.theta1_demvs<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec.theta1_demvs[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
#save(vec.theta1_demvs,file="vec_theta1_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1_demvs.Rda")
mase3 = mase(y[96:116],vec.theta1_demvs)
wmape3 = wmape(y[96:116],vec.theta1_demvs)
rmse3 = rmse(y[96:116],vec.theta1_demvs)
mae3 = mae(y[96:116],vec.theta1_demvs)

fcst_plot_4 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,fc.f = vec.theta1_demvs,timeframe = timeframe)
#save(fcst_plot_4,file="fcst_plot_4.Rda")

# TVP-reg
#--------------------------------------------------------

X = cbind(matrix(1,ncol=1,nrow=dim(X)[1]),X)

# THETA = 0.2

vec.theta02_ns<-c()
qvec.theta02_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X,activate=T)
  vec.theta02_ns[j]<-mean(first.try$fc.f)
  qvec.theta02_ns[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
#save(vec.theta02_ns,file="vec_theta02_ns.Rda")
#save(qvec.theta02_ns,file="qvec_theta02_ns.Rda")
mase(y[96:116],vec.theta02_ns)
wmape(y[96:116],vec.theta02_ns)
rmse(y[96:116],vec.theta02_ns)
mae(y[96:116],vec.theta02_ns)
lpds(y[96:116],vec.theta02_ns,qvec.theta02_ns)

fcst_plot_5 = Forecastplot(y,from=86,to=116,start_fc=96,ub=25,lb=0,
                           fc.f = vec.theta02_ns,fc.Q = qvec.theta02_ns,timeframe = timeframe)
#save(fcst_plot_5,file="fcst_plot_5.Rda")

# THETA 1

vec.theta1_ns<-c()
qvec.theta1_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec.theta1_ns[j]<-mean(second.try$fc.f)
  qvec.theta1_ns[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(vec.theta1_ns,file="vec_theta1_ns.Rda")
#save(qvec.theta1_ns,file="qvec_theta1_ns.Rda")
mase(y[96:116],vec.theta1_ns)
wmape(y[96:116],vec.theta1_ns)
rmse(y[96:116],vec.theta1_ns)
mae(y[96:116],vec.theta1_ns)
lpds(y[96:116],vec.theta1_ns,qvec.theta1_ns)

fcst_plot_6 = Forecastplot(y,from=86,to=116,start_fc=96,ub=25,lb=0,
                           fc.f = vec.theta02_ns,fc.Q = qvec.theta02_ns,
                           timeframe = timeframe)
#save(fcst_plot_6,file="fcst_plot_6.Rda")

# THETA 0.2

vec.theta02_demvs_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec.theta02_demvs_ns[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(vec.theta02_demvs_ns,file="vec_theta02_demvs_ns.Rda")
mase6 = mase(y[96:116],vec.theta02_demvs_ns)
wmape6 = wmape(y[96:116],vec.theta02_demvs_ns)
rmse6 = rmse(y[96:116],vec.theta02_demvs_ns)
mae6 = mae(y[96:116],vec.theta02_demvs_ns)

fcst_plot_7 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,
                           fc.f=vec.theta02_demvs_ns,timeframe=timeframe)
#save(fcst_plot_7,file="fcst_plot_7.Rda")


# THETA 1

vec.theta1_demvs_ns<-c()
j=1
for(i in 95:115){
  
  old_y<-y[1:i]
  old_X<-X[1:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec.theta1_demvs_ns[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
#save(vec.theta1_demvs_ns,file="vec_theta1_demvs_ns.Rda")
mase7 = mase(y[96:116],vec.theta1_demvs_ns)
wmape7 = wmape(y[96:116],vec.theta1_demvs_ns)
rmse7 = rmse(y[96:116],vec.theta1_demvs_ns)
mae7 = mae(y[96:116],vec.theta1_demvs_ns)

fcst_plot_8 = Forecastplot(y,from=86,to=116,start_fc=96,lb=0,ub=25,
                           fc.f=vec.theta1_demvs_ns,timeframe=timeframe)
#save(fcst_plot_8,file="fcst_plot_8.Rda")

#....................................................................
# Break 
#....................................................................

# from the half to the end of the sample

X = X[,2:dim(X)[2]]

# THETA = 0.2

vec2.theta02<-c()
qvec2.theta02<-c()
j=1
for(i in 183:203){
  
  old_y<-y[100:i]
  old_X<-X[100:i,]
  new_X<-X[i+1,]
  #first.try = BSTS_SeasTrend(y,S=11,U=2,sig2u=0.1,sig2d=0.1,sig2tau=0.1,N=100,n0=10,d0=10,delta=0.9,v0=0.25)
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec2.theta02[j]<-mean(first.try$fc.f)
  qvec2.theta02[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
#save(vec2.theta02,file="vec2_theta02.Rda")
#save(qvec2.theta02,file="qvec2_theta02.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta02.Rda")
mase(y[184:204],vec2.theta02)
wmape(y[184:204],vec2.theta02)
rmse(y[184:204],vec2.theta02)
mae(y[184:204],vec2.theta02)


fcst_plot_11 = Forecastplot(y,from=174,to=204,start_fc=184,lb=-30,ub=60,vec2.theta02,qvec2.theta02,timeframe)
save(fcst_plot_11,file="fcst_plot_11.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/fcst_plot_11.Rda")

# THETA 1

vec2.theta1<-c()
qvec2.theta1<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,U=2,S=11,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec2.theta1[j]<-mean(second.try$fc.f)
  qvec2.theta1[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(vec2.theta1,file="vec2_theta1.Rda")
#save(qvec2.theta1,file="qvec2_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta1.Rda")
mase(y[184:204],vec2.theta1)
wmape(y[184:204],vec2.theta1)
rmse(y[184:204],vec2.theta1)
mae(y[184:204],vec2.theta1)

fcst_plot_12 = Forecastplot(y,from=174,to=204,start_fc=184,lb=-45,ub=60,vec2.theta1,qvec2.theta1,timeframe)
save(fcst_plot_12,file="fcst_plot_12.Rda")

# THETA 0.2

vec2.theta02_demvs<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec2.theta02_demvs[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(vec2.theta02_demvs,file="vec2_theta02_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02_demvs.Rda")
mase13 = mase(y[184:204],vec2.theta02_demvs)
wmape13 = wmape(y[184:204],vec2.theta02_demvs)
rmse13 = rmse(y[184:204],vec2.theta02_demvs)
mae13 = mae(y[184:204],vec2.theta02_demvs)

fcst_plot_13 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,fc.f = vec2.theta02_demvs,timeframe = timeframe)
#save(fcst_plot_13,file="fcst_plot_13.Rda")


# THETA 1

vec2.theta1_demvs<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS.Monthly(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                           lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec2.theta1_demvs[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}
#save(vec2.theta1_demvs,file="vec2_theta1_demvs.Rda")
#load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1_demvs.Rda")
mase14 = mase(y[184:204],vec2.theta1_demvs)
wmape14 = wmape(y[184:204],vec2.theta1_demvs)
rmse14 = rmse(y[184:204],vec2.theta1_demvs)
mae14 = mae(y[184:204],vec2.theta1_demvs)

fcst_plot_14 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,fc.f = vec2.theta1_demvs,timeframe = timeframe)
#save(fcst_plot_14,file="fcst_plot_14.Rda")

# TVP-reg
#--------------------------------------------------------

X = cbind(matrix(1,ncol=1,nrow=dim(X)[1]),X)

# THETA = 0.2

vec2.theta02_ns<-c()
qvec2.theta02_ns<-c()
j=1
for(i in  183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  first.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                        phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                        lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X,activate=T)
  vec2.theta02_ns[j]<-mean(first.try$fc.f)
  qvec2.theta02_ns[j]<-mean(first.try$fc.Q)
  
  j=j+1
  
}
#save(vec2.theta02_ns,file="vec2_theta02_ns.Rda")
#save(qvec2.theta02_ns,file="qvec2_theta02_ns.Rda")
mase15 = mase(y[184:204],vec2.theta02_ns)
wmape15 = wmape(y[184:204],vec2.theta02_ns)
rmse15 = rmse(y[184:204],vec2.theta02_ns)
mae15 = mae(y[184:204],vec2.theta02_ns)


fcst_plot_15 = Forecastplot(y,from=174,to=204,start_fc=184,ub=60,lb=-30,
                            fc.f = vec2.theta02_ns,fc.Q = qvec2.theta02_ns,timeframe = timeframe)
save(fcst_plot_15,file="fcst_plot_15.Rda")

# THETA 1

vec2.theta1_ns<-c()
qvec2.theta1_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  second.try<-DSSBSTS_SV(y=old_y,X=old_X,N=100,
                         phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                         lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec2.theta1_ns[j]<-mean(second.try$fc.f)
  qvec2.theta1_ns[j]<-mean(second.try$fc.Q)
  
  j=j+1
  
}
#save(vec2.theta1_ns,file="vec2_theta1_ns.Rda")
#save(qvec2.theta1_ns,file="qvec2_theta1_ns.Rda")
mase16 = mase(y[184:204],vec2.theta1_ns)
wmape16 = wmape(y[184:204],vec2.theta1_ns)
rmse16 = rmse(y[184:204],vec2.theta1_ns)
mae16 = mae(y[184:204],vec2.theta1_ns)

fcst_plot_16 = Forecastplot(y,from=174,to=204,start_fc=184,ub=60,lb=-30,
                            fc.f = vec2.theta1_ns,fc.Q = qvec2.theta1_ns,
                            timeframe = timeframe)
save(fcst_plot_16,file="fcst_plot_16.Rda")

# THETA 0.2

vec2.theta02_demvs_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  third.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=0.2,new.X=new_X)
  vec2.theta02_demvs_ns[j]<-mean(third.try$fc.f)
  
  j=j+1
  
}
#save(vec2.theta02_demvs_ns,file="vec2_theta02_demvs_ns.Rda")
mase17 = mase(y[184:204],vec2.theta02_demvs_ns)
wmape17 = wmape(y[184:204],vec2.theta02_demvs_ns)
rmse17 = rmse(y[184:204],vec2.theta02_demvs_ns)
mae17 = mae(y[184:204],vec2.theta02_demvs_ns)

fcst_plot_17 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,
                            fc.f=vec2.theta02_demvs_ns,timeframe=timeframe)
#save(fcst_plot_17,file="fcst_plot_7.Rda")


# THETA 1

vec2.theta1_demvs_ns<-c()
j=1
for(i in 183:203){
  
  old_y<-y[90:i]
  old_X<-X[90:i,]
  new_X<-X[i+1,]
  forth.try<-DEMVS_PF(y=old_y,X=old_X,N=500,phi0=0,phi1=0.98,
                      lambda0=0.01,lambda1=0.1,THETA=1,new.X=new_X)
  vec2.theta1_demvs_ns[j]<-mean(forth.try$fc.f)
  
  j=j+1
  
}

#save(vec2.theta1_demvs_ns,file="vec2_theta1_demvs_ns.Rda")
mase18 = mase(y[184:204],vec2.theta1_demvs_ns)
wmape18 = wmape(y[184:204],vec2.theta1_demvs_ns)
rmse18 = rmse(y[184:204],vec2.theta1_demvs_ns)
mae18 = mae(y[184:204],vec2.theta1_demvs_ns)

fcst_plot_18 = Forecastplot(y,from=174,to=204,start_fc=184,lb=0,ub=25,
                            fc.f=vec2.theta1_demvs_ns,timeframe=timeframe)
#save(fcst_plot_18,file="fcst_plot_8.Rda")

#------------------------------------------------------------------
# End of forecasting
#------------------------------------------------------------------

# Break

#-------------------------------------------------------------------
# Start Tables Building
#-------------------------------------------------------------------

# Load saved results

load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta02.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta1.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta02_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/qvec2_theta1_ns.Rda")

load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta02_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec_theta1_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1_demvs.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta02_demvs_ns.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/vec2_theta1_demvs_ns.Rda")

# Generate matrix to report results via kable

list_vec = list(vec.theta02,vec.theta1,vec.theta02_ns,vec.theta1_ns,
                vec.theta02_demvs,vec.theta1_demvs,vec.theta02_demvs_ns,
                vec.theta1_demvs_ns)

list_qvec = list(qvec.theta02,qvec.theta1,qvec.theta02_ns,qvec.theta1_ns)

list_vec2 = list(vec2.theta02,vec2.theta1,vec2.theta02_ns,vec2.theta1_ns,
                 vec2.theta02_demvs,vec2.theta1_demvs,vec2.theta02_demvs_ns,
                 vec2.theta1_demvs_ns)

list_qvec2 = list(qvec2.theta02,qvec2.theta1,qvec2.theta02_ns,qvec2.theta1_ns)

MATRIX = matrix(NA,ncol = 8,nrow = 8)
n = 1
for(i in list_vec){
  
  MATRIX[n,1] = rmse(y[96:116],i)
  MATRIX[n,2] = wmape(y[96:116],i)
  MATRIX[n,3] = mase(y[96:116],i)
  n=n+1
  
}

MATRIX[1,4] = lpds(y[96:116],vec.theta02,qvec.theta02)
MATRIX[2,4] = lpds(y[96:116],vec.theta1,qvec.theta1)
MATRIX[3,4] = lpds(y[96:116],vec.theta02_ns,qvec.theta02_ns)
MATRIX[4,4] = lpds(y[96:116],vec.theta1_ns,qvec.theta1_ns)


n = 1
for(i in list_vec2){
  
  MATRIX[n,5] = rmse(y[184:204],i)
  MATRIX[n,6] = wmape(y[184:204],i)
  MATRIX[n,7] = mase(y[184:204],i)
  n=n+1
  
}

MATRIX[1,8] = lpds(y[184:204],vec2.theta02,qvec2.theta02)
MATRIX[2,8] = lpds(y[184:204],vec2.theta1,qvec2.theta1)
MATRIX[3,8] = lpds(y[184:204],vec2.theta02_ns,qvec2.theta02_ns)
MATRIX[4,8] = lpds(y[184:204],vec2.theta1_ns,qvec2.theta1_ns)

save(MATRIX,file="MATRIX.Rda")  

#----------------------------------------------------------
# End Tables building
#----------------------------------------------------------

# Break

#-------------------------------------------------------------
# Make an illustrative plot
#------------------------------------------------------------

unemp_ser = scale(LAG0$UNRATENSA)
gtrend_insurance = scale(Z_t$unemployment.insurance)
gtrend_agency = scale(Z_t$Unemployment.agency)
gtrend_depression = scale(Z_t$unemployment.depression)
gtrend_compensation = scale(Z_t$Unemployment.compensation)
gtrend_sub = scale(Z_t$subsidies)
gtrend_benefits = scale(Z_t$unemployment.benefits)
gtrend_employee = scale(Z_t$employee.benefits)
gtrend_un = scale(Z_t$unemployment)


plot_df = data.frame(timeframe,unemp_ser,gtrend_agency,gtrend_insurance,
                     gtrend_depression,gtrend_compensation,gtrend_sub,
                     gtrend_benefits,gtrend_employee,gtrend_un)

unemp_plot = ggplot(plot_df, aes(x=timeframe))+
  geom_line(aes(y=unemp_ser,col="Unemployment rate"))+
  geom_line(aes(y=gtrend_agency,col="Google searches for unemployment insurance"))+
  scale_color_manual("",
                     values=c("Unemployment rate" = "black",
                              "Google searches for unemployment insurance"= "blue"))+
  labs(x="Time",y="")+
  theme_bw()+
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
        legend.background = element_rect(fill = "white", colour = "gray30"))

save(unemp_plot,file="unemp_plot.Rda")

unemp_totalplot = ggplot(plot_df, aes(x=timeframe))+
  geom_line(aes(y=gtrend_agency,col="unemployment agency"))+
  geom_line(aes(y=gtrend_insurance,col="unemployment insurance"))+
  geom_line(aes(y=gtrend_depression,col="unemployment depression"))+
  geom_line(aes(y=gtrend_un,col="unemployment"))+
  geom_line(aes(y=unemp_ser),col="black")+
  scale_color_manual("Google searches for:",
                     values=c("unemployment insurance"= "lightgreen",
                              "unemployment agency"= "lightsalmon",
                              "unemployment depression"= "lightblue",
                              "unemployment"= "lightcoral"
                     ))+
  labs(x="Time",y="")+
  theme_bw()+
  theme(legend.direction = "vertical", legend.position = "right")

save(unemp_totalplot,file="unemp_totalplot.Rda")
