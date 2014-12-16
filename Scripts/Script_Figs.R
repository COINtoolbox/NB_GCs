# Script to prepare Figures for the paper

library(ggplot2)
library(MASS)     # glm.nb
library(COUNT)   # diagnostics
library(lmtest) # Likelihood test
library(gamlss.tr)# zero truncated model
#library(glmmADMB)  # GLMM
library(caret) # 10-fold
library(reshape) # melt
library(scales)
require(ggthemes)


GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
dim(GCS)


# Data exploration
G0<-GCS[,c("N_GC","MBH","Mdyn","Type","sig_e")]
G1<-melt(G0,id=c("N_GC","Type"))
varnames<-list(expression(log(M[BH]/M[sun])),expression(log(M[dyn]/M[sun])),expression(sigma[e]*(km/s)))
G_labeller <- function(variable,value){
  return(varnames[value])
}
ggplot(G1,aes(x=value,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_tableau()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+facet_grid(~variable, scales="free_x",labeller=G_labeller)+
  xlab("")


