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

dim(GCS[,c("MBH","upMBH","lowMBH")])
isCensored = (GCS$MBH == 0 )

N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
err_sig_e<-GCS$err_sig_e


jags.data <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = upMBH
)


model.NB.c <- "model{
  # Priors for regression coefficients
  beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
# Prior for size
size~dunif(0.001,5)
# Hyperpriors
meanx ~ dgamma(30,3)
varx ~ dgamma(2,1)
for (i in 1:N){
# MBHtrue[i]~dunif(5,12)
# MBHtrue[i]~dnorm(8,0.000001) # this would be sensible too
MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)T(5,12)
}
# Likelihood function
for (i in 1:N){
MBH[i]~dnorm(MBHtrue[i],1/errMBH[i]^2);
errorN[i]~dbin(0.5,2*errN_GC[i])
eta[i]<-beta.0+beta.1*MBHtrue[i]+exp(errorN[i]-errN_GC[i])
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


model.NB.d <- "model{
  # Priors for regression coefficients
  beta.0~dnorm(0,0.000001)
  beta.1~dnorm(0,0.000001)
  # Prior for size
  size~dunif(0.001,5)
  # Hyperpriors
  meanx ~ dgamma(30,3)
  varx ~ dgamma(2,1)
  for (i in 1:N){
    # MBHtrue[i]~dunif(5,12)
    # MBHtrue[i]~dnorm(8,0.000001) # this would be sensible too
    MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)T(5,12)
  }
  # Likelihood function
  for (i in 1:N){
    MBH[i]~dnorm(MBHtrue[i],1/errMBH[i]^2);
    errorN[i]~dbin(0.5,2*errN_GC[i])
    eta[i]<-beta.0+beta.1*MBHtrue[i]+exp(errorN[i]-errN_GC[i])
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

inits <- list(beta.0=0,beta.1=0,size=0.1)
params <- c("beta.0","beta.1","size","prediction.NB","MBHtrue")



jags.neg.c <- jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(model.NB.c),
  n.chains = 3,
  n.adapt=1000
)

jags.neg.d <- jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(model.NB.d),
  n.chains = 3,
  n.adapt=1000
)


jagssamples.nb.c <- jags.samples(jags.neg.c, params, n.iter = 1000)
jagssamples.nb.d <- jags.samples(jags.neg.d, params, n.iter = 1000)

names(jagssamples.nb.a)


summary(as.mcmc.list(jagssamples.nb.c$beta.0))
summary(as.mcmc.list(jagssamples.nb.c$beta.1))
sMBHtrue.c = summary(as.mcmc.list(jagssamples.nb.c$MBHtrue))
MBHtrue.c = sMBHtrue.c$statistics[,"Mean"]
predN.c = summary(as.mcmc.list(jagssamples.nb.c$prediction.NB))$statistics[,"Mean"]

summary(as.mcmc.list(jagssamples.nb.d$beta.0))
summary(as.mcmc.list(jagssamples.nb.d$beta.1))
sMBHtrue.d = summary(as.mcmc.list(jagssamples.nb.d$MBHtrue))
MBHtrue.d = sMBHtrue.d$statistics[,"Mean"]
predN.d = summary(as.mcmc.list(jagssamples.nb.d$prediction.NB))$statistics[,"Mean"]


points(GCS$MBH, MBHtrue.c, col="green")
points(GCS$MBH, MBHtrue.d, col="orange")
abline(a=0,b=1)

s = order(GCS$MBH)
plot(GCS$MBH[s], log(GCS$N_GC)[s])
points(GCS$MBH[s], log(predN.c)[s], col="green", type="o")
points(GCS$MBH[s], log(predN.d)[s], col="orange", type="o")

plot(MBHtrue.c, log(predN.c))
plot(MBHtrue.d, log(predN.d))






