mro_df <- read.table('random_community.size2-100.mro.table',sep = "\t",header = T)
mip_df <- read.table('random_community.size2-100.mip.table',sep = "\t",header = T)

mro_df$variable <- factor(mro_df$variable,levels = c("rand_mro","baci_a_mro","bact_mro","baci_mro","acti_mro","pseu_mro","spir_mro"))
mip_df$variable <- factor(mip_df$variable,levels = c("rand_mip","baci_a_mip","bact_mip","baci_mip","acti_mip","pseu_mip","spir_mip"))

p1 <- ggplot(mro_df)+
  geom_line(aes(pop_num,value,color=variable))+
  geom_ribbon(aes(x = pop_num,y = value, ymin = lower_ci, ymax = upper_ci, fill = variable), alpha = 0.2) +
  labs(x="Number of populations",y="MRO score")+
  scale_color_manual(values=c("#bebebe","#3ca6cc","#F8f16f","#5bb288","#ab84b6","#ff4f42","#BC587D"))+
  scale_fill_manual(values=c("#bebebe","#3ca6cc","#F8f16f","#5bb288","#ab84b6","#ff4f42","#BC587D"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  expand_limits(x=0,y=0.5)+
  scale_x_continuous(limits = c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits = c(0.4,0.7),breaks=seq(0.4,0.7,0.1))  

p2 <- ggplot(mip_df)+
  geom_line(aes(pop_num,value,color=variable))+
  geom_ribbon(aes(x = pop_num,y = value, ymin = lower_ci, ymax = upper_ci, fill = variable), alpha = 0.2) +
  labs(x="Number of populations",y="MIP score")+
  scale_color_manual(values=c("#bebebe","#3ca6cc","#F8f16f","#5bb288","#ab84b6","#ff4f42","#BC587D"))+
  scale_fill_manual(values=c("#bebebe","#3ca6cc","#F8f16f","#5bb288","#ab84b6","#ff4f42","#BC587D"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  expand_limits(x=0,y=120)+
  scale_x_continuous(limits = c(0,100),breaks = seq(0,100,25))+
  scale_y_continuous(limits = c(0,120),breaks=seq(0,120,40)) 

library(gridExtra)
plots <- list(p1,p2)
grid.arrange(grobs = plots,ncol=2)
