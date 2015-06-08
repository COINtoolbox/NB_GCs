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
dim(GCS)
N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
err_sig_e<-GCS$err_sig_e
N = nrow(GCS)


# Data list for JAGS
jags.data <- list(
  N_GC = GCS$N_GC,
  MBH = GCS$MBH,
  N = nrow(GCS)
)

# Poisson model

model.pois<-"model{
#Priors for regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)

# Likelihood function
for (i in 1:N){
eta[i]<-beta.0+beta.1*MBH[i]
log(mu[i])<-max(-20,min(20,eta[i]))# Ensures that large beta values do not cause numerical problems. 
N_GC[i]~dpois(mu[i])
#
# Prediction
prediction.pois[i]~dpois(mu[i])
# Discrepancy measures
YNew[i] ~ dpois(mu[i])
expY[i] <- mu[i]
varY[i] <- mu[i]
PRes[i] <-(N_GC[i] - expY[i])/sqrt(varY[i])
PResNew[i] <-(YNew[i] - expY[i])/sqrt(varY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)
}
Fit<-sum(D[1:N])
New<-sum(DNew[1:N])
}"

#inits<-list(beta.0=coefficients(glm.pois)[1],beta.1=coefficients(glm.pois)[2])
inits<-list(beta.0=0,beta.1=0)

params<-c("beta.0","beta.1","prediction.pois","PRes")

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
#summary(as.mcmc.list(jagssamples$beta.0))
#summary(as.mcmc.list(jagssamples$beta.1))

N_low<-GCS$N_GC-N_err

N_low[N_low<0]<-0
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}
CairoPDF("..//Figures/JAGS_pois.pdf",height=8,width=9)
ggplot(pred.pois2,aes(x=MBH,y=NGC))+
  geom_ribbon(aes(ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
  geom_ribbon(aes(ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray") +
  geom_ribbon(aes(ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray") +
  geom_point(aes(colour=Type,shape=Type),size=3.25)+
  geom_errorbar(guide="none",aes(colour=Type,ymin=N_low,ymax=NGC+N_err),alpha=0.7,width=0.05)+
  geom_errorbarh(guide="none",aes(colour=Type,xmin=MBH-GCS$lowMBH,
                                  xmax=MBH+upMBH),alpha=0.7,height=0.05)+
  geom_line(aes(x=MBH,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_y_continuous(trans = 'asinh',breaks=c(0,10,100,1000,10000,100000),labels=c("0",expression(10^1),expression(10^2),
                                                                                   expression(10^3),expression(10^4),expression(10^5)))+
  scale_colour_gdocs()+
  scale_shape_manual(values=c(19,2,8))+
  theme_hc()+
  theme_hc()+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[BH]/M['\u0298']))+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()


require(scales)
Pres<-summary(as.mcmc.list(jagssamples$PRes),quantiles=0.5)$quantiles
Dipersion = sum(Pres^2)/(N-2)# beta.0, beta.1, 2 parameters

