#Poisson and NB regression using JAGS by Rafael S. de Souza, Bart Buelens, Ewan Cameron

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


# Function to allow parse labels in facet_wrap
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}
give.n <- function(x){
  
  return(c(y = 0.5, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
################

# Read data

GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
dim(GCS)
N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
err_sig_e<-GCS$err_sig_e




######## NB with errors ########################################################

jags.data3 <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = upMBH
)

model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,5)

# Hyperpriors

meanx ~ dgamma(30,3)
varx ~ dgamma(2,1)

for (i in 1:N){

#MBHtrue[i]~dunif(5,12)

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
#prediction.NB[i]~dnegbin(p[i],size)
}
}"
inits3 <- list(beta.0=0,beta.1=0,size=0.1)
params3 <- c("beta.0","beta.1","size","prediction.NB","MBHtrue")

jags.neg3 <- jags.model(
  data = jags.data3, 
  inits = inits3, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=1000
)

update(jags.neg3, 10000)

jagssamples.nb3 <- jags.samples(jags.neg3, params3, n.iter = 5000)
codasamples.nb3 <- coda.samples(jags.neg3, params3, n.iter = 5000)


ggs(as.mcmc.list(jagssamples.nb3) ,family=c("beta"))

summary(as.mcmc.list(jagssamples.nb3$beta.0))
summary(as.mcmc.list(jagssamples.nb3$beta.1))
summary(as.mcmc.list(jagssamples.nb3$size))

MBHtrue<-summary(as.mcmc.list(jagssamples.nb3$MBHtrue),quantiles=0.5)
pred.NBerr<-summary(as.mcmc.list(jagssamples.nb3$prediction.NB),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MBHtrue=MBHtrue$quantiles,MBH=GCS$MBH,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])

S.NB1<-ggs(codasamples.nb3 ,family=c("beta"))
S.NB2<-ggs(codasamples.nb3,family=c("size"))
S.NB<-rbind(S.NB1,S.NB2,deparse.level=2)
S.NB$Parameter<-revalue(S.NB$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1]),
                                          "size"="k"))

ggs_density(S.NB)+
  scale_colour_economist(guide="none")+
  theme_hc()+
  scale_fill_economist()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))+xlab("Parameter  value")+ylab("Density")



CairoPDF("JAGS_NB.pdf",height=8,width=9)
ggplot(pred.NB2err,aes(x=MBH,y=NGC))+
  geom_ribbon(aes(x=MBHtrue,y=mean,ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(aes(x=MBHtrue,y=mean,ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(aes(x=MBHtrue,y=mean,ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_err,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
                                  xmax=MBH+upMBH),alpha=0.7)+
  geom_line(aes(x=MBHtrue,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'log10',breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[BH]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()


#Not accounting for errors ########################################################
jags.data<-list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  N = nrow(GCS)
)

# Poisson model

model.pois<-"model{
# Priors for regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
# Likelihood function
for (i in 1:N){
eta[i]<-beta.0+beta.1*MBH[i]
log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems. 
N_GC[i]~dpois(mu[i])
# Prediction
prediction.pois[i]~dpois(mu[i])
}

}"

#inits<-list(beta.0=coefficients(glm.pois)[1],beta.1=coefficients(glm.pois)[2])
inits<-list(beta.0=0,beta.1=0)

params<-c("beta.0","beta.1","prediction.pois")

jags.pois<-jags.model(
  data = jags.data, 
  inits = inits, 
  textConnection(model.pois),
  n.chains = 3,
  n.adapt=1000
  
)
update(jags.pois, 20000)
posterior.pois <- coda.samples(jags.pois, params, n.iter = 50000)
jagssamples <- jags.samples(jags.pois, params, n.iter = 50000)
pred.pois<-summary(as.mcmc.list(jagssamples$prediction.pois),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.pois2<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MBH=GCS$MBH,mean=pred.pois$quantiles[,4],lwr1=pred.pois$quantiles[,3],lwr2=pred.pois$quantiles[,2],lwr3=pred.pois$quantiles[,1],upr1=pred.pois$quantiles[,5],upr2=pred.pois$quantiles[,6],upr3=pred.pois$quantiles[,7])

# Posterior means of beta for comparison with ML estimates
summary(as.mcmc.list(jagssamples$beta.0))
summary(as.mcmc.list(jagssamples$beta.1))

CairoPDF("JAGS_pois.pdf",height=8,width=9)
ggplot(pred.pois2,aes(x=MBH,y=NGC))+
  geom_ribbon(aes(ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(aes(ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(aes(ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_err,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
                                  xmax=MBH+upMBH),alpha=0.7)+
  geom_line(aes(x=MBH,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'log10',breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[BH]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()


S.pois<-ggs(posterior.pois,family="beta")
S.pois$Parameter<-revalue(S.pois$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1])))


# Plot 

g0<-ggs_traceplot(S.pois)+
  scale_colour_economist()+
  #  theme_pander(base_size = 20,nomargin = F)+
  theme_hc()+scale_alpha_manual(values=c(0.3,0.3,0.3))+
  geom_line(alpha=0.5)+scale_linetype_manual(values=c("solid","dotted","dashed"))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab("Parameter value")+
  xlab("Iteration")+theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
                          legend.position="none",plot.title = element_text(hjust=0.5),
                          axis.title.y=element_text(vjust=0.75),
                          axis.title.x=element_text(vjust=-0.25),
                          text = element_text(size=25))+
  facet_grid(Parameter~.,labeller=label_parsed,scales = "free")

CairoPDF("chain_poisson.pdf",height=10,width=8)
g0 
dev.off()


g1<-ggs_density(S.pois)+
  scale_colour_economist(guide="none")+
  theme_hc()+
  scale_fill_economist()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))+xlab("Parameter  value")+ylab("Density")
CairoPDF("posterior_poisson.pdf",height=10,width=8)
facet_wrap_labeller(g1,labels=c(expression(beta[0]),expression(beta[1])))
dev.off()


# Negative Binomial version


model.NB<-"model{
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
# Prediction
prediction.NB[i]~dnegbin(p[i],size)
}

}"

inits2<-list(beta.0=0,beta.1=0,size=0.1)
params2<-c("beta.0","beta.1","size","prediction.NB")

jags.neg<-jags.model(
  data = jags.data, 
  inits = inits2, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=1000
  
)
update(jags.neg, 20000)
posterior.NB <- coda.samples(jags.neg, params2, n.iter = 50000)
jagssamples.nb <- jags.samples(jags.neg, params2, n.iter = 50000)
pred.NB<-summary(as.mcmc.list(jagssamples.nb$prediction.NB),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))

pred.NB2<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MBH=GCS$MBH,mean=pred.NB$quantiles[,4],lwr1=pred.NB$quantiles[,3],lwr2=pred.NB$quantiles[,2],lwr3=pred.NB$quantiles[,1],upr1=pred.NB$quantiles[,5],upr2=pred.NB$quantiles[,6],upr3=pred.NB$quantiles[,7])

# Posterior means of beta for comparison with ML estimates
summary(as.mcmc.list(jagssamples.nb$beta.0))
summary(as.mcmc.list(jagssamples.nb$beta.1))

CairoPDF("JAGS_NB.pdf",height=8,width=9)
ggplot(pred.NB2,aes(x=MBH,y=NGC))+
  geom_ribbon(aes(ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(aes(ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(aes(ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_err,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
                                  xmax=MBH+upMBH),alpha=0.7)+
  geom_line(aes(x=MBH,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'log10',breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[BH]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()


gelman.diag(posterior.NB)
S.NB1<-ggs(posterior.NB,family=c("beta"))
S.NB2<-ggs(posterior.NB,family=c("size"))
S.NB<-rbind(S.NB1,S.NB2,deparse.level=2)
S.NB$Parameter<-revalue(S.NB$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1]),
                                          "size"="k"))


p0<-ggs_traceplot(S.NB)+
  scale_colour_economist()+
  #  theme_pander(base_size = 20,nomargin = F)+
  theme_hc()+scale_alpha_manual(values=c(0.3,0.3,0.3))+
  geom_line(alpha=0.5)+scale_linetype_manual(values=c("solid","dotted","dashed"))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab("Parameter value")+
  xlab("Iteration")+theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
                          legend.position="none",plot.title = element_text(hjust=0.5),
                          axis.title.y=element_text(vjust=0.75),
                          axis.title.x=element_text(vjust=-0.25),
                          text = element_text(size=25))+
  facet_grid(Parameter~.,labeller=label_parsed,scales = "free")

CairoPDF("chain_NB",height=10,width=8)
p0
#facet_wrap_labeller(p0,labels=c(expression(beta[0]),expression(beta[1]))) 
dev.off()

p1<-ggs_density(S.NB)+
  scale_colour_economist(guide="none")+
  theme_hc()+
  scale_fill_economist()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))+xlab("Parameter  value")+ylab("Density")
CairoPDF("posterior_NB.pdf",height=10,width=8)
facet_wrap_labeller(p1,labels=c(expression(beta[0]),expression(beta[1]),"k"))
dev.off()

ggs_ppmean(ggs(posterior.NB,family=c("prediction")),  outcome=GCS$N_GC)


ggs_ppsd(ggs(posterior.NB,family=c("prediction")),  outcome=GCS$N_GC)
