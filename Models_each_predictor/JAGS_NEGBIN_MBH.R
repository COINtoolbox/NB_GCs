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
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
#err_sig_e<-GCS$err_sig_e


######## NB with errors ########################################################
MBHx = seq(from = 0.95 * min(GCS$MBH), 
           to = 1.05 * max(GCS$MBH), 
           length.out = 100)

jags.data3 <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  errMBH = upMBH,
  MBHx = MBHx,
  M = 100
)

model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,5)

# Hyperpriors

#meanx ~ dgamma(30,3)
#varx ~ dgamma(2,1)
meanx ~ dgamma(85,10)
varx ~ dgamma(2,1)

for (i in 1:N){

#MBHtrue[i]~dunif(5,12)

# MBHtrue[i]~dnorm(8,0.000001) # this would be sensible too
MBHtrue[i] ~ dgamma(meanx^2/varx,meanx/varx)T(6,11)
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
FitNew<-sum(DNew[1:N])
# Prediction for new data
for (j in 1:M){
  etax[j]<-beta.0+beta.1*MBHx[j]
  log(mux[j])<-max(-20,min(20,etax[j]))
  px[j]<-size/(size+mux[j])
  prediction.NBx[j]~dnegbin(px[j],size)
}
}"
inits3 <- list(beta.0=0,beta.1=0,size=0.1)
params3 <- c("beta.0","beta.1","size","prediction.NB","MBHtrue","Fit","FitNew","prediction.NBx")

jags.neg3 <- jags.model(
  data = jags.data3, 
  inits = inits3, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=1000
)

update(jags.neg3, 10000)

jagssamples.nb3 <- jags.samples(jags.neg3, params3, n.iter = 50000)
codasamples.nb3 <- coda.samples(jags.neg3, params3, n.iter = 50000)
dicsamples.nb3 <- dic.samples(jags.neg3, params3, n.iter = 50000,type="pD")



summary(as.mcmc.list(jagssamples.nb3$beta.0))
summary(as.mcmc.list(jagssamples.nb3$beta.1))
summary(as.mcmc.list(jagssamples.nb3$size))

MBHtrue<-summary(as.mcmc.list(jagssamples.nb3$MBHtrue),quantiles=0.5)
pred.NBerr<-summary(as.mcmc.list(jagssamples.nb3$prediction.NB),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MBHtrue=MBHtrue$quantiles,MBH=GCS$MBH,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])
pred.NBerrx<-summary(as.mcmc.list(jagssamples.nb3$prediction.NBx),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2errx<-data.frame(MBHx=MBHx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])





#CairoPDF("..//Figures/JAGS_NBx.pdf",height=8,width=9)
#CairoFonts(regular = 'Calibri:style=Regular')
CairoPDF("..//Figures/JAGS_NBx.pdf",height=8,width=9)
ggplot(pred.NB2err,aes(x=MBH,y=NGC))+
  geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(data=pred.NB2errx,aes(x=MBHx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_err,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
                                  xmax=MBH+upMBH),alpha=0.7)+
  geom_line(data=pred.NB2errx,aes(x=MBHx,y=mean),colour="gray25",linetype="dashed",size=1.2)+
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


CairoPDF("..//Figures/JAGS_NB_MBH.pdf",height=8,width=9)
#pdf("..//Figures/JAGS_NB_MBH.pdf",height=8,width=9)
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



S.NB1<-ggs(codasamples.nb3 ,family=c("beta"))
S.NB2<-ggs(codasamples.nb3,family=c("size"))
#S.NB3<-ggs(codasamples.nb3,family=c("Fit"))


# Diagnostics
Pred<-ggs(codasamples.nb3,family=c("FitNew"))[,"value"]
Obs<-ggs(codasamples.nb3,family=c("Fit"))[1:75000,"value"]
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
















