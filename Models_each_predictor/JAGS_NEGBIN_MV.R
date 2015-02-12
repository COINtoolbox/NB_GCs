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





# Script starts here 


# Read data

GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
#dim(GCS)
N_err<-GCS$N_GC_err
err_MV_T<-GCS$err_MV_T


######## NB with errors ########################################################

jags.data3 <- list(
  N_GC = GCS$N_GC,
  MV_T = GCS$MV_T,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  err_MV_T = err_MV_T
  # values for gamma hyperpriors  
 )


model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,5)

#

for (i in 1:N){

MV_T_true[i]~dunif(-25,-15)
}

# Likelihood function

for (i in 1:N){

MV_T[i]~dnorm(MV_T_true[i],1/err_MV_T[i]^2);

errorN[i]~dbin(0.5,2*errN_GC[i])

eta[i]<-beta.0+beta.1*MV_T_true[i]+exp(errorN[i]-errN_GC[i])

log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems.

p[i]<-size/(size+mu[i])

N_GC[i]~dnegbin(p[i],size)

# Prediction
etaTrue[i]<-beta.0+beta.1*MV_T_true[i]
    log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
    pTrue[i]<-size/(size+muTrue[i])
prediction.NB[i]~dnegbin(pTrue[i],size)
#prediction.NB[i]~dnegbin(p[i],size)

# Discrepancy measures
YNew[i] ~ dnegbin(p[i],size)
expY[i] <- mu[i]
varY[i] <- mu[i] + pow(mu[i],2) / size
PRes[i] <-(N_GC[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

}
Fit<-sum(D[1:N])
NewFit<-sum(DNew[1:N])
}"
inits3 <- list(beta.0=0,beta.1=0,size=0.1)
params3 <- c("beta.0","beta.1","size","prediction.NB","MV_T_true","Fit","NewFit")

jags.neg3 <- jags.model(
  data = jags.data3, 
  inits = inits3, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=1000
)

update(jags.neg3, 10000)

jagssamples.nb3 <- jags.samples(jags.neg3, params3, n.iter = 10000)
codasamples.nb3 <- coda.samples(jags.neg3, params3, n.iter = 10000)
dicsamples.nb3 <- dic.samples(jags.neg3, params3, n.iter = 10000,type="pD")



summary(as.mcmc.list(jagssamples.nb3$beta.0))
summary(as.mcmc.list(jagssamples.nb3$beta.1))
summary(as.mcmc.list(jagssamples.nb3$size))

MV_T_true<-summary(as.mcmc.list(jagssamples.nb3$MV_T_true),quantiles=0.5)
pred.NBerr<-summary(as.mcmc.list(jagssamples.nb3$prediction.NB),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MV_T_true=MV_T_true$quantiles,MV_T=GCS$MV_T,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])

S.NB1<-ggs(codasamples.nb3 ,family=c("beta"))
S.NB2<-ggs(codasamples.nb3,family=c("size"))
#S.NB3<-ggs(codasamples.nb3,family=c("Fit"))


# Diagnostics
Pred<-ggs(codasamples.nb3,family=c("NewFit"))[,"value"]
Obs<-ggs(codasamples.nb3,family=c("Fit"))[1:30000,"value"]
sqrt(mean((Pred-Obs)^2))




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



pdf("..//Figures/JAGS_NB_M_V.pdf",height=8,width=9)
ggplot(pred.NB2err,aes(x=MV_T,y=NGC))+
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_err,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MV_T-GCS$err_MV_T,
                                  xmax=MV_T+err_MV_T),alpha=0.7)+
  geom_line(aes(x=MV_T_true,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'log10',breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(M[V]))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()





















