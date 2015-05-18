#NB regression using JAGS by Rafael S. de Souza, Bart Buelens, Ewan Cameron and Joseph Hilbe

#  Required libraries
library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(plyr)
library(MASS)
library(scales)
require(runjags)


# Script starts here 


# Read data

GCS = read.csv(file="GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
N = nrow(GCS)

######## NB with errors ########################################################
MBHx = seq(from = 0.95 * min(GCS$MBH), 
           to = 1.05 * max(GCS$MBH), 
           length.out = 500)

jags.data <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = upMBH,
  MBHx = MBHx,
  M = 500
)

model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,10)

# Hyperpriors

meanx ~ dgamma(0.01,0.01)
varx ~ dgamma(0.01,0.01)

for (i in 1:N){


MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)
}

# Likelihood function

for (i in 1:N){

MBH[i]~dnorm(MBHtrue[i],1/errMBH[i]^2);

errorN[i]~dbin(0.5,2*errN_GC[i])

eta[i]<-beta.0+beta.1*MBHtrue[i]

log(mu[i])<-log(exp(eta[i])+errorN[i]-errN_GC[i])

#Using mixture of Poisson and Gamma

N_GC[i]~dpois(g[i])
g[i]~dgamma(size,rateParm[i])
rateParm[i]<-size/mu[i]

# Discrepancy measures
YNew[i] ~ dpois(g[i])
expY[i] <- mu[i]
varY[i] <- mu[i] + pow(mu[i],2) / size
PRes[i] <-(N_GC[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

}
Fit<-sum(D[1:N])
New<-sum(DNew[1:N])

# Prediction for new data
for (j in 1:M){
etax[j]<-beta.0+beta.1*MBHx[j]
log(mux[j])<-max(-20,min(20,etax[j]))
prediction.NBx[j]~dpois(gx[j])
gx[j]~dgamma(size,rateParmx[j])
rateParmx[j]<-size/mux[j]
}
}"
inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
params <- c("beta.0","beta.1","size","PRes","prediction.NB","MBHtrue","Fit","New","prediction.NBx")

library(parallel)
cl <- makeCluster(3)
jags.neg <- run.jags(method="rjparallel", method.options=list(cl=cl),
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model.NB,
                     n.chains = 3,
                     adapt=2000,
                     monitor=c(params),
                     burnin=20000,
                     sample=50000,
                     summarise=FALSE,
                     plots=FALSE
)
jagssamples.nb <- as.mcmc.list(jags.neg )
summary<-extend.jags(jags.neg,drop.monitor=c("PRes","prediction.NB","MBHtrue","Fit","New","prediction.NBx"), summarise=TRUE)
print(summary)


