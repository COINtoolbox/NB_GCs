# compare different error models for negbin

library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(plyr)
library(MASS)
library(scales)

GCS = read.csv(file="..//Dataset//GCs_full.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(MBH)) # 1 removed
#UpM<-min(GCS$MBH)
UpM<-GCS$MBH+GCS$upMBH
#UpM<-3
# Censoring information
isCensored = (GCS$MBH == 0 )
GCS$MBH[isCensored] = NA
#censorLimitVec = rep(UpM, length(GCS$MBH) )



treshMat<-matrix(rep(c(quantile(UpM[isCensored]),max(UpM)),length(GCS$MBH)),nrow=length(GCS$MBH),ncol=6,byrow=TRUE)
censorLimitVec=UpM[isCensored]


xbin<-censorLimitVec
xbin[xbin<=treshMat[1,1]]<-0
xbin[xbin<=treshMat[1,2] & xbin > treshMat[1,1]]<-1
xbin[xbin<=treshMat[1,3]  & xbin > treshMat[1,2] ]<-2
xbin[xbin<=treshMat[1,4]  & xbin > treshMat[1,3]]<-3
xbin[xbin<=treshMat[1,5]  & xbin > treshMat[1,4]]<-4
xbin[xbin<=treshMat[1,6]  & xbin > treshMat[1,5]]<-5

MBHbin<-rep(5,length(GCS$MBH) )
MBHbin[isCensored]<-xbin

# Initialize 
xinit=rep( NA , length(GCS$MBH) )
for(i in 1:length(GCS$MBH)){
  if (is.na(GCS$MBH[i])){
    if (MBHbin[i]==0){
      xinit[i]=treshMat[i,1]-1
    } else if (MBHbin[i]==4){
      xinit[i]=treshMat[i,ncol(treshMat)-1]+1
  } else {
    xinit[i]= treshMat[i,3]
  }
}
}
 



N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
err_sig_e<-GCS$err_sig_e


jags.data <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = upMBH,
  
# isCensored = as.numeric(isCensored), 
 MBHbin = MBHbin, 
 treshMat = treshMat
 )


model.NB.c <- "model{
  # Priors for regression coefficients
  beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
# Prior for size
size~dunif(0.001,5)
# Hyperpriors
#meanx ~ dgamma(30,3)
#varx ~ dgamma(2,1)

meanx ~ dgamma(0.5,0.5)
varx ~ dgamma(0.5,0.5)
for (i in 1:N){

MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)

}
# Likelihood function
for (i in 1:N){
MBHbin[i] ~ dinterval(MBHtrue[i],  treshMat[i,])
MBH[i]~dnorm(MBHtrue[i],1/errMBH[i]^2);

errorN[i]~dbin(0.5,2*errN_GC[i])
eta[i]<-beta.0+beta.1*MBHtrue[i]
log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems.
p[i]<-size/(size+mu[i])
N_GC[i]~dnegbin(p[i],size)

# Prediction
etaTrue[i]<-beta.0+beta.1*MBHtrue[i]
log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
pTrue[i]<-size/(size+muTrue[i])
prediction.NB[i]~dnegbin(pTrue[i],size)
}
}"

inits <- list(beta.0=0,beta.1=0,size=0.1,MBHtrue=xinit)
params <- c("beta.0","beta.1","size","prediction.NB","MBHtrue")


jags.neg.c <- jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(model.NB.c),
  n.chains = 3,
  n.adapt=1000
)











jagssamples.nb.c <- jags.samples(jags.neg.c, params, n.iter = 1000)



summary(as.mcmc.list(jagssamples.nb.c$beta.0))
summary(as.mcmc.list(jagssamples.nb.c$beta.1))
sMBHtrue.c = summary(as.mcmc.list(jagssamples.nb.c$MBHtrue))
MBHtrue.c = sMBHtrue.c$statistics[,"Mean"]
predN.c = summary(as.mcmc.list(jagssamples.nb.c$prediction.NB))$statistics[,"Mean"]




points(GCS$MBH, MBHtrue.c, col="green")
abline(a=0,b=1)

s = order(GCS$MBH)
plot(GCS$MBH[s], log(GCS$N_GC)[s])
points(GCS$MBH[s], log(predN.c)[s], col="green", type="o")


plot(MBHtrue.c, log(predN.c))







