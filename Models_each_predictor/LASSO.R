#  R script  GLMM
#  Copyright (C) 2014  Rafael S. de Souza, Bart Buelens, Ewan Cameron
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

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

GCS = read.csv(file="..//Dataset//GCs_full.csv",header=TRUE,dec=".",sep="")
Full_type<-read.table(file="..//Dataset//fulltype.dat",header=TRUE)
GCS$alltype<-Full_type$fulltype
GCS$alltype<-GCS$alltype
GCS = subset(GCS, !is.na(MV_T)) 
#dim(GCS)
N_err<-GCS$N_GC_err
err_MV_T<-GCS$err_MV_T
N = nrow(GCS)
#type<-match(GCS$alltype,unique(GCS$alltype))

type<-match(Full_type$fulltype,unique(Full_type$fulltype))
Ntype<-length(unique(GCS$alltype))






######## NB with errors ########################################################
MV_Tx = seq(from = 1.05 * min(GCS$MV_T), 
            to = 0.95 * max(GCS$MV_T), 
            length.out = 500)

jags.data <- list(
  N_GC = GCS$N_GC,
  MV_T = GCS$MV_T,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  err_MV_T = err_MV_T,
  MV_Tx = MV_Tx,
  M = 500,
  type=type,
  Ntype=Ntype
)


model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)
beta.2~dnorm(0,0.000001)

#tau.R~dgamma(0.01,0.01)
tau.R<-pow(sdBeta,-1)
#sdBeta ~ dunif(0,20)
sdBeta ~ dgamma(0.01,0.01)


# Prior for size

size~dunif(0.001,5)

#

for (i in 1:N){

MV_T_true[i]~dunif(-26,-10)
}

for (j in 1:Ntype){

ranef[j]~ddexp(0,tau.R)
}


# Likelihood function

for (i in 1:N){

MV_T[i]~dnorm(MV_T_true[i],1/err_MV_T[i]^2);

errorN[i]~dbin(0.5,2*errN_GC[i])

eta[i]<-beta.0+beta.1*MV_T_true[i]+ranef[type[i]]

log(mu[i])<-log(exp(eta[i])+errorN[i]-errN_GC[i])

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
New<-sum(DNew[1:N])

# Prediction for new data
for (j in 1:M){
etax[j]<-beta.0+beta.1*MV_Tx[j]
log(mux[j])<-max(-20,min(20,etax[j]))
px[j]<-size/(size+mux[j])
prediction.NBx[j]~dnegbin(px[j],size)
}
}"
inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))

params <- c("beta.0","beta.1","size","ranef","PRes","MV_T_true","Fit","New","prediction.NBx","tau.R")

library(parallel)
cl <- makeCluster(3)
jags.neg <- run.jags(method="rjparallel", method.options=list(cl=cl),
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model.NB,
                     n.chains = 3,
                     adapt=1000,
                     monitor=c(params),
                     burnin=20000,
                     sample=30000,
                     summarise=FALSE,
                     thin=2,
                     plots=FALSE
)

jagssamples.nb <- as.mcmc.list(jags.neg )
summary<-extend.jags(jags.neg,drop.monitor=c("PRes","MV_T_true","Fit","New","prediction.NBx","tau.R"), summarise=TRUE)


MV_T_true<-summary(as.mcmc.list(jags.neg,vars="MV_T_true"),quantiles=0.5)
#pred.NBerr<-summary(as.mcmc.list(jagssamples.nb, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
#pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,MV_T_true=MV_T_true$quantiles,MV_T=GCS$MV_T,mean=pred.NBerr$statistics[,1],lwr1=pred.NBerr$quantiles[,3],lwr2=pred.NBerr$quantiles[,2],lwr3=pred.NBerr$quantiles[,1],upr1=pred.NBerr$quantiles[,5],upr2=pred.NBerr$quantiles[,6],upr3=pred.NBerr$quantiles[,7])
pred.NBerrx<-summary(as.mcmc.list(jags.neg, vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2errx<-data.frame(MV_Tx=MV_Tx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])



N_low<-GCS$N_GC-N_err

N_low[N_low<0]<-0


asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

cairo_pdf("..//Figures/M_Vx_random.pdf",height=8,width=9)
ggplot(GCS,aes(x=MV_T,y=N_GC))+
  geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray",method = "loess") +
  geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray",method = "loess") +
  geom_ribbon(data=pred.NB2errx,aes(x=MV_Tx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray",method = "loess") +
  geom_point(aes(colour=Type,shape=Type),size=3.25,alpha=0.8)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=N_low,ymax=N_GC+N_err),alpha=0.7,width=0.05)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MV_T-GCS$err_MV_T,
                                  xmax=MV_T+err_MV_T),alpha=0.7,height=0.05)+
  geom_line(data=pred.NB2errx,aes(x=MV_Tx,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),
                                                                                   expression(10^3),expression(10^4),expression(10^5)))+
  
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8,10))+scale_x_reverse()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(M[V]))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                               axis.title.y=element_text(vjust=0.75),
                               axis.title.x=element_text(vjust=-0.25),
                               text = element_text(size=25))
dev.off()



L.radon.intercepts <- data.frame(
  Parameter=paste("ranef[", seq(1:69), "]", sep=""),
  Label=levels(Full_type$fulltype))
head(L.radon.intercepts)

S.full <- ggs(jagssamples.nb,par_labels=L.radon.intercepts,family=c("ranef"))
library(RColorBrewer)
blues_fun <- colorRampPalette(brewer.pal(9,"Blues")[4:9])
blues=blues_fun(69)


pdf("..//Figures/random_LASSO.pdf",height=14,width=9)
ggs_caterpillar(S.full)+geom_vline(aes(yintercept=0),color="gray80",size=1,linetype="dashed")+
  theme_hc()+ 
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
         axis.title.y=element_text(size=25,vjust=0.75),
         axis.title.x=element_text(size=25,vjust=-0.25),axis.text.x =element_text(size=25),
         text = element_text(size=17))+aes(color=Parameter)+
  scale_color_manual(guide="none",values = blues)+ylab("Galaxy Type")+
  xlab(expression(paste(zeta[i]," Highest Posterior Density"," ")))
dev.off()



pdf("..//Figures/JAGS_NB_M_V.pdf",height=8,width=9)
ggplot(pred.NB2err,aes(x=MV_T,y=NGC))+
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(aes(x=MV_T_true,y=mean,ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=NGC-N_low,ymax=NGC+N_err),alpha=0.7)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MV_T-GCS$err_MV_T,
                                  xmax=MV_T+err_MV_T),alpha=0.7)+
  geom_line(aes(x=MV_T_true,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),
                                                                                   expression(10^3),expression(10^4),expression(10^5)))+
  
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8,10))+scale_x_reverse()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(M[V]))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                               axis.title.y=element_text(vjust=0.75),
                               axis.title.x=element_text(vjust=-0.25),
                               text = element_text(size=25))
dev.off()



# Diagnostics






S.NB1<-ggs(jagssamples.nb ,family=c("beta"))
S.NB2<-ggs(jagssamples.nb,family=c("size"))

S.NB<-rbind(S.NB1,S.NB2,deparse.level=2)
S.NB$Parameter<-revalue(S.NB$Parameter, c("beta.0"=expression(beta[0]), "beta.1"=expression(beta[1]),
                                          "size"="k"))


g1<-ggs_density(S.NB)+
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
CairoPDF("..//Figures/posterior_MV_full.pdf",height=10,width=8)
facet_wrap_labeller(g1,labels=c(expression(beta[0]),expression(beta[1]),"k"))
dev.off()





# Dispersion parameter

require(scales)
Pres<-summary(as.mcmc.list(jags.neg, vars="PRes"),quantiles=0.5)$quantiles
Dipersion = sum(Pres^2)/(N-72)# beta.0, beta.1 and k, 3 parameters + 69 random intercepts



# Model comparison 
Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
sqrt(mean((Pred-Obs)^2))


dicsamples.nb <- dic.samples(jags.neg, params, n.iter = 50000,type="pD")



# Plot residuals vc galaxy type
clus_data<-data.frame(Pres=Pres,type=GCS$alltype)
p <- ggplot(clus_data, aes(x=type, y=Pres),group=type)+ xlab("Galaxy Type") +
  ylab("Pearson Residuals")





pdf("..//Figures/Pres_random_LASSO.pdf",height=6,width=14)
p + stat_boxplot(colour="gray",geom ='errorbar')+geom_boxplot(aes(group=type,colour=type,fill=type),outlier.shape = 19,colour="gray",fatten=2,size=1,outlier.size=2,outlier.colour = "gray",notchwidth = 0.35,notch=F,data=clus_data)+
  theme_hc()+
  scale_fill_manual(guide="none",values = blues)+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(angle=-90,size=12.5),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))
dev.off()









