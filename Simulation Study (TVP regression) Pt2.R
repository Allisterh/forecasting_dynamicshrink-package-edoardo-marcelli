# Simulation Study (TVP regression) Pt2 with Precision sampler
#---------------------------------------------

prova = DSS_PRECISIONSAMPLER(y=y,X=X,N=1800,burn=200,phi0=0,phi1=0.98,gamma0=0.5,v0=0.25,
                             lambda0=0.01,lambda1=0.1,THETA=0.2)

#save(prova,file="chan_1000.Rda")
load("C:/Users/edoar/Dropbox/elementi salvati R/chan_1000.Rda")

time = c(1:length(y))
j = 4
m = prova$beta[,j]
C.m = prova$var.beta[,j]
trueb.m = true_b[,j]
dataframe = data.frame(m,C.m,time,trueb.m)

ggplot(dataframe, aes(x=time))+
  geom_line(aes(y=trueb.m))+
  geom_line(aes(y=m),col="orange")+
  coord_cartesian(ylim = c(-4, +4))+
  geom_ribbon(aes(ymin=m-1.96*sqrt(C.m),ymax=m+1.96*sqrt(C.m)),alpha=0.16,fill="orange")+
  labs(x="Time",
       y="")+
  theme_bw()

save(plot_chan_1,file="plot_chan_1.Rda")
save(plot_chan_2,file="plot_chan_2.Rda")
save(plot_chan_3,file="plot_chan_3.Rda")
save(plot_chan_4,file="plot_chan_4.Rda")
save(plot_chan_5,file="plot_chan_5.Rda")
save(plot_chan_6,file="plot_chan_6.Rda")

v = exp(prova$h)
dataframe = data.frame(v,time)

plot_chan_v <- ggplot(dataframe, aes(x=time))+
  geom_line(aes(y=v))+
  #geom_ribbon(aes(ymin=v-1.96*sqrt(var.v),ymax=v+1.96*sqrt(var.v)),alpha=0.26,fill="grey")+
  # geom_ribbon(aes(ymin=v-1.282*sqrt(var.v),ymax=v+1.282*sqrt(var.v)),alpha=0.16,fill="grey")+
  #geom_ribbon(aes(ymin=v-2.576*sqrt(var.v),ymax=v+2.576*sqrt(var.v)),alpha=0.16,fill="grey")+
  coord_cartesian(ylim = c(0, 1))+
  labs(x="Time",
       y="")+
  theme_bw()

#save(plot_chan_v,file="plot_chan_v.Rda")

j=6
ind_chan = prova$ind[,j]
ind_df = MCMC.out(df.50,200)$ind[,j]
ind_sv = MCMC.out(tvp.reg.sv.50,200)$ind[,j]
ind_true = true_ind[,j]

dataframe = data.frame(ind_chan,ind_df,ind_true,ind_sv,time)

plot_ind_6 = ggplot(dataframe, aes(x=time))+
  geom_line(aes(y=ind_chan),col="red")+
  geom_line(aes(y=ind_df),col="blue")+
  geom_line(aes(y=ind_sv),col="green")+
  geom_line(aes(y=ind_true), linetype="dotted")+
  coord_cartesian(ylim = c(0, +1))+
  labs(x="Time",
       y="")+
  theme_bw()

save(plot_ind_1,file="plot_ind_1.Rda")
save(plot_ind_2,file="plot_ind_2.Rda")
save(plot_ind_3,file="plot_ind_3.Rda")
save(plot_ind_4,file="plot_ind_4.Rda")
save(plot_ind_5,file="plot_ind_5.Rda")
save(plot_ind_6,file="plot_ind_6.Rda")


estimated.theta = prova$beta
estimated.gammino = prova$ind

SMM_MAT = matrix(NA,ncol=6,nrow=1)
SMM_MAT[,1] = SSE(true_b,estimated.theta)
SMM_MAT[,2] = HmD(true_ind,estimated.gammino)
SMM_MAT[,3] = False.Positive(true_ind,estimated.gammino)
SMM_MAT[,4] = False.Negative(true_ind,estimated.gammino)
SMM_MAT[,5] = False.Active(true_ind,estimated.gammino)
SMM_MAT[,6] = False.Non_Active(true_ind,estimated.gammino)

#save(SMM_MAT,file="SMM_MAT.Rda")

SMM_MAT_2 = matrix(NA,ncol=4,nrow=1)
SMM_MAT_2[,1] = SSE(true_b[,1:4],estimated.theta[,1:4])
SMM_MAT_2[,2] = HmD(true_ind[,1:4],estimated.gammino[,1:4])
SMM_MAT_2[,3] = SSE(true_b[,5:50],estimated.theta[,5:50])
SMM_MAT_2[,4] = HmD(true_ind[,5:50],estimated.gammino[,5:50])

#save(SMM_MAT_2,file="SMM_MAT_2.Rda")

