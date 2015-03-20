# Model comparison
require(ggplot2)
require(ggthemes)
require(reshape)
require(Cairo)
Models<-as.factor(c("MBH","Mdyn","MV","sigma"))
DIC<-c(628.9,687.7,685.8,1015)
RMS<-c(24.6085,24.98447,22.4822,28.25993)
Dispersion<-c(1.034658,1.243808,0.8932604,1.201396)
comp<-data.frame(Models,DIC,RMS,Dispersion)


c <- ggplot(comp, aes(x=Models,y=DIC, fill=Models,colour=Models))
CairoPDF("..//Figures/DIC.pdf",height=8,width=9)
c + geom_bar(stat = "identity") + 
#  facet_wrap(~test,scales="free")+
  coord_flip()+
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
  scale_x_discrete(breaks=c("MBH","Mdyn","MV","sigma"),labels=c(expression(M[BH]),expression(M[dyn]),expression(M[V]),expression(sigma)))+
  xlab("")
dev.off()



c2 <- ggplot(comp, aes(x=Models,y=RMS, fill=Models,colour=Models))
CairoPDF("..//Figures/RMS.pdf",height=8,width=9)
c2 + geom_bar(stat = "identity") + 
  #  facet_wrap(~test,scales="free")+
  coord_flip()+
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
  scale_x_discrete(breaks=c("MBH","Mdyn","MV","sigma"),labels=c(expression(M[BH]),expression(M[dyn]),expression(M[V]),expression(sigma)))+
  xlab("")
dev.off()




c3 <- ggplot(comp, aes(x=Models,y=Dispersion-1, fill=Models,colour=Models))
CairoPDF("..//Figures/Dispersion_stat.pdf",height=8,width=9)
c3 + geom_bar(origin = 1,stat = "identity")   +
  #  facet_wrap(~test,scales="free")+
  coord_flip()+
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
  scale_x_discrete(breaks=c("MBH","Mdyn","MV","sigma"),labels=c(expression(M[BH]),expression(M[dyn]),expression(M[V]),expression(sigma)))+
  ylab("Dispersion Statistics")+geom_hline(yintercept = 1,linetype="dotted",colour="red",size=1.5)+xlab("")+
  scale_y_continuous(limits=c(-0.3, 0.3), breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),labels=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)+1)
dev.off()
