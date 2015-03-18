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

GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
N_err<-GCS$N_GC_err
N = nrow(GCS)


######## NB with errors ########################################################
Mdynx = seq(from = 0.95 * min(GCS$Mdyn), 
           to = 1.05 * max(GCS$Mdyn), 
           length.out = 500)

jags.data <- list(
  N_GC = GCS$N_GC,
  Mdyn = GCS$Mdyn,
  errN_GC = GCS$N_GC_err,
  N = nrow(GCS),
  Mdynx = Mdynx,
  M = 500
  )

model.NB <- "model{

# Priors for regression coefficients

beta.0~dnorm(0,0.000001)

beta.1~dnorm(0,0.000001)

# Prior for size

size~dunif(0.001,5)


# Likelihood function

for (i in 1:N){


errorN[i]~dbin(0.5,2*errN_GC[i])

eta[i]<-beta.0+beta.1*Mdyn[i]

#log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems.
log(mu[i])<-log(exp(eta[i])+errorN[i]-errN_GC[i])
p[i]<-size/(size+mu[i])

N_GC[i]~dnegbin(p[i],size)

# Prediction
etaTrue[i]<-beta.0+beta.1*Mdyn[i]
  log(muTrue[i])<-max(-20,min(20,etaTrue[i]))
    pTrue[i]<-size/(size+muTrue[i])
prediction.NB[i]~dnegbin(pTrue[i],size)


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
  etax[j]<-beta.0+beta.1*Mdynx[j]
  log(mux[j])<-max(-20,min(20,etax[j]))
  px[j]<-size/(size+mux[j])
  prediction.NBx[j]~dnegbin(px[j],size)
}
}"
#inits <- list(beta.0=0,beta.1=0,size=0.1)
inits1 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits2 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
inits3 <- list(beta.0=rnorm(1,0,0.1),beta.1=rnorm(1,0,0.1),size=runif(1,0.1,5))
params <- c("beta.0","beta.1","size","PRes","prediction.NB","Fit","New","prediction.NBx")


library(parallel)
cl <- makeCluster(3)
jags.neg <- run.jags(method="rjparallel", method.options=list(cl=cl),
                     data = jags.data, 
                     inits = list(inits1,inits2,inits3),
                     model=model.NB,
                     n.chains = 3,
                     adapt=2000,
                     monitor=c(params),
                     burnin=15000,
                     sample=25000,
                     summarise=FALSE,
                     plots=FALSE
)

jagssamples.nb <- as.mcmc.list(jags.neg )


#jags.neg <- jags.model(
#  data = jags.data, 
#  inits = inits, 
#  textConnection(model.NB),
#  n.chains = 3,
#  n.adapt=1000
#)

#update(jags.neg, 10000)

#jagssamples.nb <- jags.samples(jags.neg, params, n.iter = 50000)





#pred.NBerr<-summary(as.mcmc.list(jags.neg, vars="prediction.NB"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
#pred.NB2err<-data.frame(Type=GCS$Type,NGC=GCS$N_GC,Mdyn=GCS$Mdyn,mean=pred.NBerr$statistics[1:N,1],lwr1=pred.NBerr$quantiles[1:N,3],lwr2=pred.NBerr$quantiles[1:N,2],lwr3=pred.NBerr$quantiles[1:N,1],upr1=pred.NBerr$quantiles[1:N,5],upr2=pred.NBerr$quantiles[1:N,6],upr3=pred.NBerr$quantiles[1:N,7])
pred.NBerrx<-summary(as.mcmc.list(jags.neg,vars="prediction.NBx"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
pred.NB2errx<-data.frame(Mdynx=Mdynx,mean=pred.NBerrx$statistics[,1],lwr1=pred.NBerrx$quantiles[,3],lwr2=pred.NBerrx$quantiles[,2],lwr3=pred.NBerrx$quantiles[,1],upr1=pred.NBerrx$quantiles[,5],upr2=pred.NBerrx$quantiles[,6],upr3=pred.NBerrx$quantiles[,7])


asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

#N_low<-asinh(pred.NB2err$NGC-N_err)
N_low<-GCS$N_GC-N_err

N_low[N_low<0]<-0

#CairoPDF("..//Figures/Mdynx.pdf",height=8,width=9)
CairoPDF("..//Figures/Mdynx.pdf",height=8,width=9)
ggplot(GCS,aes(x=Mdyn,y=N_GC))+
  geom_ribbon(data=pred.NB2errx,aes(x=Mdynx,y=mean,ymin=lwr1, ymax=upr1), alpha=0.3, fill="gray") +
  geom_ribbon(data=pred.NB2errx,aes(x=Mdynx,y=mean,ymin=lwr2, ymax=upr2), alpha=0.2, fill="gray") +
  geom_ribbon(data=pred.NB2errx,aes(x=Mdynx,y=mean,ymin=lwr3, ymax=upr3), alpha=0.1, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=N_low,ymax=N_GC+N_err),alpha=0.7,width=0.05)+
  #  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
  #                                  xmax=MBH+upMBH),alpha=0.7)+
  geom_line(data=pred.NB2errx,aes(x=Mdynx,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),
                                                                                   expression(10^3),expression(10^4),expression(10^5)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[dyn]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                 axis.title.y=element_text(vjust=0.75),
                                                 axis.title.x=element_text(vjust=-0.25),
                                                 text = element_text(size=25))
dev.off()




# Diagnostics

#Density 

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
CairoPDF("..//Figures/posterior_Mdyn.pdf",height=10,width=8)
facet_wrap_labeller(g1,labels=c(expression(beta[0]),expression(beta[1]),"k"))
dev.off()

# Mixing of chains 

g0<-ggs_traceplot(S.NB)+
  scale_colour_economist(guide="none")+
  theme_hc()+
  scale_fill_economist()+
  #  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))+
  ylab("Parameter value")+
  xlab("Iteration")+
  facet_grid(Parameter~.,labeller=label_parsed,scales = "free")

CairoPDF("chain_Mdyn.pdf",height=10,width=8)
g0 
dev.off()



# Model comparison 

Pred<-ggs(jagssamples.nb,family=c("New"))[,"value"]
Obs<-ggs(jagssamples.nb,family=c("Fit"))[,"value"]
sqrt(mean((Pred-Obs)^2))


# Dispersion parameter

require(scales)
Pres<-summary(as.mcmc.list(jags.neg, vars="PRes"),quantiles=0.5)$quantiles
Dipersion = sum(Pres^2)/(N-3)# beta.0, beta.1 and k, 3 parameters



# Plot residuals vc galaxy type
clus_data<-data.frame(Pres=Pres,Mdyn=GCS$Mdyn,type=GCS$Type)
p <- ggplot(clus_data, aes(x=type, y=Pres),group=type)+ xlab("Galaxy Type") +
  ylab("Pearson Residuals")

pdf("Pres_Mdyn.pdf",height=5.5,width=9)
p + stat_boxplot(colour="gray",geom ='errorbar')+geom_boxplot(aes(group=type,colour=type,fill=type),outlier.shape = 19,colour="gray",fatten=2,size=1,outlier.size=2,outlier.colour = "gray",notchwidth = 0.35,notch=F,data=clus_data)+
  theme_hc()+
  scale_fill_economist()+
  theme(strip.background = element_rect(fill="gray95"),plot.background = element_rect(fill = 'white', colour = 'white'),
        legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25))
dev.off()




# DIC 

jags.DIC <- jags.model(
  data = jags.data, 
  inits = inits1, 
  textConnection(model.NB),
  n.chains = 3,
  n.adapt=2000
)

update(jags.DIC , 10000)
dicsamples.nb <- dic.samples(jags.DIC, params, n.iter = 25000,type="pD")
















