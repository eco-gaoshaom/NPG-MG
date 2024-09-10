library(dplyr)
library(ggplot2)
library(ggbreak)

pri_data <- read.delim("site.rand_mro_mip_div",header = T,row.names = 1)
pri_data <- pri_data %>%  mutate(mro_new = if_else(mro > rad_mro + rad_mro_ci, mro, -mro)) 
pri_data <- pri_data %>%  mutate(mip_new = if_else(mip > rad_mip + rad_mip_ci, mip, -mip)) 
pri_data$intertype <- ifelse(pri_data$mro_new > 0 & pri_data$mip_new > 0, "Type 1",ifelse(pri_data$mro_new < 0 & pri_data$mip_new > 0, "Type 2",ifelse(pri_data$mro_new < 0 & pri_data$mip_new < 0, "Type 3","Type 4")))
pri_data$diet_habitat <- paste(pri_data$diet_new,pri_data$habitat_new,sep = "_")
pri_data$diet_habitat <- factor(pri_data$diet_habitat,
  levels = c("non-folivore_Wild","folivore_Wild","non-folivore_Captive","folivore_Captive"))
pri_data$group <- factor(pri_data$group,levels = c("Ape","Lemur","Monkey"))

pri_data_t1 <- pri_data[pri_data$intertype == "Type 1",]
pri_data_t2 <- pri_data[pri_data$intertype == "Type 2",]
pri_data_t3 <- pri_data[pri_data$intertype == "Type 3",]
pri_data_t4 <- pri_data[pri_data$intertype == "Type 4",]

#X Range:-0.55 to -0.45 and 0.52 to 0.62
range(pri_data_t1$mro_new)
range(pri_data_t4$mro_new) 
range(pri_data_t2$mro_new)
range(pri_data_t3$mro_new)

#Y Range:-150 to -30 and 40 to 160
range(pri_data_t1$mip_new)
range(pri_data_t2$mip_new) 
range(pri_data_t3$mip_new)
range(pri_data_t4$mip_new)

#geom_point for type1
p2 <- ggplot(pri_data_t1,aes(mro_new,mip_new))+
  geom_point(size=2.5,aes(shape=diet_habitat,color = group),alpha = 0.75)+
  scale_shape_manual(values = c(1,2,16,17))+
  scale_color_manual(values=c("#3ca6cc","#ff4f42"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),
        legend.text = element_text(size=4),
        legend.position = 'top')+
  labs(x="Competition (MRO score)",y="Cooperation (MIP score)")+
  scale_x_continuous(limits = c(0.52,0.62),breaks=seq(0.52,0.62,0.05))+
  scale_y_continuous(limits = c(40,160),breaks=seq(40,160,60))+
  coord_cartesian(xlim= c(0.52,0.62),ylim = c(40,160))+
  expand_limits(x=0.52,y=40)


p1 <- ggplot(pri_data_t2,aes(mro_new,mip_new))+
  geom_point(size=2.5,aes(shape=diet_habitat,color = group),alpha = 0.75)+
  scale_shape_manual(values = c(1,2,16,17))+
  scale_color_manual(values=c("#3ca6cc","#ab84b6","#ff4f42"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),
        legend.text = element_text(size=4),
        legend.position = 'top')+
  labs(x="Competition (MRO score)",y="Cooperation (MIP score)")+
  scale_x_continuous(limits = c(-0.55,-0.45),breaks=seq(-0.55,-0.45,0.05))+
  scale_y_continuous(limits = c(40,160),breaks=seq(40,160,60))+
  coord_cartesian(xlim= c(-0.55,-0.45),ylim = c(40,160))+
  expand_limits(x=-0.55,y=40)

p3 <- ggplot(pri_data_t3,aes(mro_new,mip_new))+
  geom_point(size=2.5,aes(shape=diet_habitat,color = group),alpha = 0.75)+
  scale_shape_manual(values = c(1,2,16,17))+
  scale_color_manual(values=c("#3ca6cc","#ab84b6","#ff4f42"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),
        legend.text = element_text(size=4),
        legend.position = 'top')+
  labs(x="Competition (MRO score)",y="Cooperation (MIP score)")+
  scale_x_continuous(limits = c(-0.55,-0.45),breaks=seq(-0.55,-0.45,0.05))+
  scale_y_continuous(limits = c(-150,-30),breaks=seq(-150,-30,60))+
  coord_cartesian(xlim= c(-0.55,-0.45),ylim = c(-150,-30))+
  expand_limits(x=-0.55,y=-150)


p4 <- ggplot(pri_data_t4,aes(mro_new,mip_new))+
  geom_point(size=2.5,aes(shape=diet_habitat,color = group),alpha = 0.75)+
  scale_shape_manual(values = c(1,2,16,17))+
  scale_color_manual(values=c("#3ca6cc","#ab84b6","#ff4f42"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),
        legend.text = element_text(size=4),
        legend.position = 'top')+
  labs(x="Competition (MRO score)",y="Cooperation (MIP score)")+
  scale_x_continuous(limits = c(0.52,0.62),breaks=seq(0.52,0.62,0.05))+
  scale_y_continuous(limits = c(-150,-30),breaks=seq(-150,-30,60))+
  coord_cartesian(xlim= c(0.52,0.62),ylim = c(-150,-30))+
  expand_limits(x=0.52,y=-150)

library(gridExtra)
plots <- list(p1,p2,p3,p4)
grid.arrange(grobs = plots,ncol=2)
