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
require(Cairo)
require(plyr)


GCS = read.csv(file="..//Dataset//GCs.csv",header=TRUE,dec=".",sep="")
GCS = subset(GCS, !is.na(Mdyn)) # 1 removed
dim(GCS)


# Data exploration
G0<-GCS[,c("N_GC","Type","MBH","Mdyn","sig_e")]
G1<-melt(G0,id=c("N_GC","Type"))
varnames<-list(expression(log~M[BH]/M['\u0298']),expression(log~M[dyn]/M['\u0298']),
               expression(sigma[e]*(km/s)),expression(R[e]*(kpc)))
G_labeller <- function(variable,value){
  return(varnames[value])
}
N_err<-GCS$N_GC_err
lowMBH<-GCS$lowMBH
upMBH<-GCS$upMBH
err_sig_e<-GCS$err_sig_e


CairoPDF("MBH.pdf",height=8,width=9,family="Symbol")
ggplot(GCS,aes(x=MBH,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.8)+
  geom_errorbarh(aes(xmin=MBH-GCS$lowMBH,
                     xmax=MBH+upMBH),alpha=0.8)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+scale_shape_stata()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[BH]/M['\u0298']))+theme(plot.title = element_text(hjust=0.5),
                 axis.title.y=element_text(vjust=0.75),
                 axis.title.x=element_text(vjust=-0.25),
                 text = element_text(size=25))
dev.off()

CairoPDF("Mdyn.pdf",height=8,width=9)
ggplot(GCS,aes(x=Mdyn,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.8)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_gdocs()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(log~M[dyn]/M['\u0298']))+theme(plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()

CairoPDF("sig_e.pdf",height=8,width=9,family="Symbol")
ggplot(GCS,aes(x=sig_e,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.7)+
  geom_errorbarh(aes(xmin=sig_e-err_sig_e,
                     xmax=sig_e+err_sig_e),alpha=0.7)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_tableau()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(sigma~(km/s)))+theme(plot.title = element_text(hjust=0.5),
                                                axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()

CairoPDF("Re.pdf",height=8,width=9,family="Symbol")
ggplot(GCS,aes(x=Re,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.7)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_tableau()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(R[e]~(kpc)))+theme(plot.title = element_text(hjust=0.5),
                                          axis.title.y=element_text(vjust=0.75),
                                          axis.title.x=element_text(vjust=-0.25),
                                          text = element_text(size=25))
dev.off()


CairoPDF("MV.pdf",height=8,width=9,family="Symbol")
ggplot(GCS,aes(x=MV_T,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.7)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_tableau()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(M[V]))+theme(plot.title = element_text(hjust=0.5),
                                     axis.title.y=element_text(vjust=0.75),
                                     axis.title.x=element_text(vjust=-0.25),
                                     text = element_text(size=25))+
  scale_x_reverse()
dev.off()

CairoPDF("MK.pdf",height=8,width=9,family="Symbol")
ggplot(GCS,aes(x=MK,y=N_GC,colour=Type,shape=Type))+geom_point(size=3)+
  geom_errorbar(guide="none",aes(ymin=N_GC-N_err,ymax=N_GC+N_err),alpha=0.7)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_colour_tableau()+scale_shape_tremmel()+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  ylab(expression(N[GC]))+
  xlab(expression(M[K]))+theme(plot.title = element_text(hjust=0.5),
                               axis.title.y=element_text(vjust=0.75),
                               axis.title.x=element_text(vjust=-0.25),
                               text = element_text(size=25))+
  scale_x_reverse()
dev.off()

