#we firstly plot boxplot for mro/mip across different diet/habitat/primate groups
library(ggplot2)
library(data.table)
library(ggsignif)
pri_data <- read.delim("site.rand_mro_mip_div",header = T,row.names = 1)
pri_data$diet_new <- factor(pri_data$diet_new,levels = c("folivore","non-folivore"))
pri_data$habitat_new <- factor(pri_data$habitat_new,levels = c("Wild","Captive"))
pri_data$group <- factor(pri_data$group,levels = c("Ape","Lemur","Monkey"))

ggplot(pri_data,aes(diet_new,mro))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = diet_new),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = diet_new),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#90bb97","#9b87be"))+
  labs(x = "Diet group",y="Competition (MRO score)")+
  geom_signif(comparisons = list(c("folivore","non-folivore")),
              map_signif_level = F, 
              textsize = 4, 
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.45,0.65),breaks = seq(0.45,0.65,0.05))+
  coord_cartesian(ylim= c(0.45,0.65))+
  expand_limits(y=0.45)  

ggplot(pri_data,aes(habitat_new,mro))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = habitat_new),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = habitat_new),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  
  scale_color_manual(values=c("#5bb288","#905797"))+
  labs(x = "Habitat group",y="Competition (MRO score)")+
  geom_signif(comparisons = list(c("Wild","Captive")),
              map_signif_level = F, 
              textsize = 4, 
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.45,0.65),breaks = seq(0.45,0.65,0.05))+
  coord_cartesian(ylim= c(0.45,0.65))+
  expand_limits(y=0.45) 

ggplot(pri_data,aes(group,mro))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.35,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#3ca6cc","#ab84b6","#ff4f42"))+
  labs(x = "Phylogenetic group",y="Competition (MRO score)")+
  geom_signif(comparisons = list(c("Ape","Lemur"),c("Lemur","Monkey"),c("Ape","Monkey")),
              map_signif_level = F, 
              textsize = 4, 
              step_increase = 0.1,
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.45,0.65),breaks = seq(0.45,0.65,0.05))+
  coord_cartesian(ylim= c(0.45,0.65))+
  expand_limits(y=0.45)


ggplot(pri_data,aes(diet_new,mip))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = diet_new),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = diet_new),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#90bb97","#9b87be"))+
  labs(x = "Diet group",y="Cooperation (MIP score)")+
  geom_signif(comparisons = list(c("folivore","non-folivore")),
              map_signif_level = F, 
              textsize = 4, 
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0,180),breaks = seq(0,180,45))+
  coord_cartesian(ylim= c(0,180))+
  expand_limits(y=0)  

ggplot(pri_data,aes(habitat_new,mip))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = habitat_new),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = habitat_new),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#5bb288","#905797"))+
  labs(x = "Habitat group",y="Cooperation (MIP score)")+
  geom_signif(comparisons = list(c("Wild","Captive")),
              map_signif_level = F, 
              textsize = 4, 
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0,180),breaks = seq(0,180,45))+
  coord_cartesian(ylim= c(0,180))+
  expand_limits(y=0)  

ggplot(pri_data,aes(group,mip))+
  geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.35,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#3ca6cc","#ab84b6","#ff4f42"))+
  labs(x = "Phylogenetic group",y="Cooperation (MIP score)")+
  geom_signif(comparisons = list(c("Ape","Lemur"),c("Lemur","Monkey"),c("Ape","Monkey")),
              map_signif_level = F, 
              textsize = 4, 
              step_increase = 0.05,
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0,180),breaks = seq(0,180,45))+
  coord_cartesian(ylim= c(0,180))+
  expand_limits(y=0)  
