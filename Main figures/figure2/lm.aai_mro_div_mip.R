library(ggplot2)
pri_data <- read.delim("site.rand_mro_mip_div",header = T,check.names = F,row.names = 1)

summary(lm(pri_data$mro~pri_data$aai))
summary(lm(pri_data$mro~pri_data$div))
summary(lm(pri_data$mip~pri_data$aai))
summary(lm(pri_data$mip~pri_data$div))

p1 <- ggplot(pri_data)+
  geom_point(aes(div,mro),size=1,color="#c8c8c8")+
  geom_smooth(aes(div,mro),method='lm',fill=NA,color = "#5bb288")+
  labs(x="Richness",y="Competition (MRO)")+
  #scale_color_manual(values=c(rgb(194,220,80,max=255),rgb(171,132,182,max=255)))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=0,y=0.45)+
  scale_x_continuous(limits = c(0,900),breaks = seq(0,900,450))+
  scale_y_continuous(limits = c(0.45,0.65),breaks=seq(0.45,0.65,0.1))  

p2 <- ggplot(pri_data)+
  geom_point(aes(aai,mro),size=1,color="#c8c8c8")+
  geom_smooth(aes(aai,mro),method='lm',fill=NA,color = "#5bb288")+
  labs(x="Mean genomic relatedness",y="Competition (MRO)")+
  #scale_color_manual(values=c(rgb(194,220,80,max=255),rgb(171,132,182,max=255)))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=35,y=0.45)+
  scale_x_continuous(limits = c(37,43),breaks = seq(37,43,3))+
  scale_y_continuous(limits = c(0.45,0.65),breaks=seq(0.45,0.65,0.1))  

p3 <- ggplot(pri_data)+
  geom_point(aes(div,mip),size=1,color="#c8c8c8")+
  geom_smooth(aes(div,mip),method='lm',fill=NA,color = "#905797")+
  labs(x="Richness",y="Cooperation (MIP)")+
  #scale_color_manual(values=c(rgb(194,220,80,max=255),rgb(171,132,182,max=255)))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=0,y=0.085)+
  scale_x_continuous(limits = c(0,900),breaks = seq(0,900,450))+
  scale_y_continuous(limits = c(0,160),breaks = seq(0,160,80)) 

p4 <- ggplot(pri_data)+
  geom_point(aes(aai,mip),size=1,color="#c8c8c8")+
  geom_smooth(aes(aai,mip),method='lm',fill=NA,color = "#905797")+
  labs(x="Mean genomic relatedness",y="Cooperation (MIP)")+
  #scale_color_manual(values=c(rgb(194,220,80,max=255),rgb(171,132,182,max=255)))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=37,y=40)+
  scale_x_continuous(limits = c(37,43),breaks = seq(37,43,3))+
  scale_y_continuous(limits = c(0,160),breaks = seq(0,160,80)) 

library(gridExtra)
plots <- list(p1,p2,p3,p4)
grid.arrange(grobs = plots,ncol=2)
