# GLM Jags

GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
dim(GCS)

# run ML GLM to get initial values for jags
glm.pois<-glm(N_GC ~ MBH, GCS, family="poisson")
glm.neg<-glm.nb(N_GC ~ MBH, GCS)

jags.data<-list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  N = nrow(GCS)
  )

# Poisson version
k.bugs<-"model{
# Priors for regression coefficients

beta.0~dnorm(0,000001)
beta.1~dnorm(0,0.000001)
# Likelihood function
for (i in 1:N){
eta[i]<-beta.0+beta.1*MBH[i]
log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems. 
N_GC[i]~dpois(mu[i])
}

}"

inits<-list(beta.0=coefficients(glm.pois)[1],beta.1=coefficients(glm.pois)[2])
params<-c("beta.0","beta.1")

jags.m<-jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(k.bugs),
  n.chains = 3,
  n.adapt=1000
 
 )
update(jags.m, 10000)
samps <- coda.samples(jags.m, params, n.iter = 10000)
plot(samps)
gelman.diag(samps)


# Negative Binomial version

# run ML GLM to get initial values for jags

k2.bugs<-"model{
# Priors for regression coefficients
beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
# Prior for size 
size~dunif(0.001,5)

# Likelihood function
for (i in 1:N){
eta[i]<-beta.0+beta.1*MBH[i]
log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems. 
p[i]<-size/(size+mu[i])
N_GC[i]~dnegbin(p[i],size)
}

}"

inits2<-list(beta.0=coefficients(glm.neg)[1],beta.1=coefficients(glm.neg)[2],size=0.1)
params<-c("beta.0","beta.1","size")

jags.neg<-jags.model(
  data = jags.data, 
  inits = inits2, 
  textConnection(k2.bugs),
  n.chains = 3,
  n.adapt=1000
  
)
update(jags.neg, 10000)
samps2 <- coda.samples(jags.neg, params, n.iter = 50000)
plot(samps2)
gelman.diag(samps2)



