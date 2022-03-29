#------------------------------------------------------------------------#
#                                                                        #
#       Accuracy measures                                                #
#                                                                        #
#------------------------------------------------------------------------#

mae = function(y,yhat){mean(abs(y - yhat))}
me = function(y,yhat){mean(y - yhat)}
mpe = function(y,yhat){mean((y-yhat)/y)}
mape = function(y, yhat){mean(abs((y - yhat)/y))}
rmse =  function(y, yhat){sqrt(mean(((y - yhat)^2)))}
wmape = function(y, yhat) sum(abs(y-yhat))/sum(abs(y))
mase = function(y,yhat){mae(y,yhat)/mae(y[-1],head(y,-1))}
SSE = function(y,yhat){sum((y-yhat)^2)}
lpds = function(y,meanhat,varhat){
  sumvec = c()
  for(i in 1:length(y)){
    sumvec[i] = dnorm(y[i],mean=meanhat[i],sd=sqrt(varhat[i]),log=T)
  }
  sumofsumvec = sum(sumvec)
  sumofsumvec
}

HmD = function(true,estimate){
  estimate[estimate>=0.5]=1
  estimate[estimate<0.5]=0
  sum1 = c()
  for(j in 1:(dim(estimate)[2])){
    binaryvec = abs(estimate[,j] - true[,j])
    sum1[j] = sum(binaryvec)
  }
  sum2 = sum(sum1)
  sum2
}

False.Positive = function(true,estimate){
  estimate[estimate>0.5]=1
  estimate[estimate<0.5]=0
  detect.1 = matrix(0,dim(estimate)[1],dim(estimate)[2])
  detect.1[estimate == 1 & true==0]=1
  sum(detect.1)
}

False.Negative = function(true,estimate){
  estimate[estimate>0.5]=1
  estimate[estimate<0.5]=0
  detect.0 = matrix(0,dim(estimate)[1],dim(estimate)[2])
  detect.0[estimate == 0 & true==1]=1
  sum(detect.0)
}

False.Active = function(true,estimate){
  estimate[estimate>0.5]=1
  estimate[estimate<0.5]=0
  active.var = c()
  active.var[colSums(estimate)>0 & colSums(true)==0]=1
  active.var[is.na(active.var)]=0
  sum(active.var)
}

False.Non_Active = function(true,estimate){
  estimate[estimate>0.5]=1
  estimate[estimate<0.5]=0
  active.var = c()
  active.var[colSums(estimate)==0 & colSums(true)>0]=1
  active.var[is.na(active.var)]=0
  sum(active.var)
}

lpds_multi = function(y,fc_m,fc_v){
  require(mvtnorm)
  list_matv=list()
  for(i in 1:dim(y)[1]){
    v = matrix(fc_v[,i])
    mat=matrix(NA,9,9)
    mat[lower.tri(mat,diag=T)]=v
    mat[upper.tri(mat)]=mat[lower.tri(mat)]
    list_matv[[i]]=mat
  }
  
  sumvec = c()
  for(i in 1:dim(y)[1]){
    sumvec[i] = dmvnorm(y[i,],mean=fc_m[,i],sigma=list_matv[[i]],log=T)
  }
  sumofsumvec = sum(sumvec)
  sumofsumvec
}
