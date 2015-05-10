library(LaplacesDemon)
require(ggplot2)
require(ggthemes)




#Plot Probability Functions
x <- seq(from=-2, to=2, by=0.01)
Lap1<-data.frame(x,dlap=dlaplace(x,0,0.25))
Lap2<-data.frame(x,dlap=dlaplace(x,0,0.5))
Gauss1<-data.frame(x,Dnorm=dnorm(x,0,0.25))
Gauss2<-data.frame(x,Dnorm=dnorm(x,0,0.5))

cairo_pdf("..//Figures/Laplace_prior.pdf",height=8,width=9)
ggplot(Lap1,aes(x=x,y=dlap))+
  geom_line(color="blue3")+
  geom_line(data=Gauss,aes(x=x,y=Dnorm),color="blue3",linetype="dashed")+
  geom_line(data=Gauss2,aes(x=x,y=Dnorm),color="cyan3",linetype="dashed")+
  geom_line(data=Lap2,aes(x=x,y=dlap),color="cyan3")+
   theme_hc()+
  theme( axis.title.y=element_text(vjust=0.75),
         axis.title.x=element_text(vjust=-0.25),
         text = element_text(size=25))+ylab("Density distribution")
dev.off()
  


