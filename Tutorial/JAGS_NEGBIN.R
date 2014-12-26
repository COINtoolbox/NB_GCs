library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)

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

beta.0~dnorm(0,0.000001)
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
  n.chains = 4,
  n.adapt=1000
 
 )
update(jags.m, 10000)
samps <- coda.samples(jags.m, params, n.iter = 10000)
summary(samps)
ggs(samps)
S1<-ggs(samps)
S1$Parameter<-revalue(S1$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1])))

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
  n.chains = 4,
  n.adapt=1000
  
)
update(jags.neg, 10000)
samps2 <- jags.samples(jags.neg, params, n.iter = 10000)
plot(samps2)
gelman.diag(samps2)
S2<-ggs(samps2)


library(plyr)
S2$Parameter<-revalue(S2$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1]),
              "size"="k"))


# Plot 
CairoPDF("chain1.pdf")
ggs_traceplot(S1)+
  scale_colour_gdocs()+theme_pander(base_size = 15,nomargin = F)+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab("Value")+
  xlab("Iteration")+theme(plot.background = element_rect(fill = 'white', colour = 'white'),
                          legend.position="none",plot.title = element_text(hjust=0.5),
                          axis.title.y=element_text(vjust=0.75),
                          axis.title.x=element_text(vjust=-0.25),
                          text = element_text(size=25))+
  facet_grid(Parameter~.,labeller=label_parsed,scales = "free")
dev.off()

CairoPNG("chain2.png")
ggs_traceplot(S2)+
scale_colour_gdocs()+theme_pander(base_size = 15,nomargin = F)+
#  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab("Value")+
  xlab("Iteration")+theme(plot.background = element_rect(fill = 'white', colour = 'white'),
                          legend.position="none",plot.title = element_text(hjust=0.5),
                                       axis.title.y=element_text(vjust=0.75),
                                       axis.title.x=element_text(vjust=-0.25),
                                       text = element_text(size=25))+
  facet_grid(Parameter~.,labeller=label_parsed,scales = "free")
dev.off()

plot(samps)
gelman.diag(samps)
str(samps)

