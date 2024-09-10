library(dplyr)
#we calculated the mro/mip variation by deleting each bin in a community
site_bin_mro_mip_var <- read.table('site_bin_mro_mip_var',header = F)
del_var <- site_bin_mro_mip_var[,-1]
colnames(del_var) <- c("sample","bin","mro_var","mip_var")

pri_data <- read.delim("site.rand_mro_mip_div",header = T,row.names = 1)
pri_data <- pri_data %>%  mutate(mro_new = if_else(mro > rad_mro + rad_mro_ci, mro, -mro)) 
pri_data <- pri_data %>%  mutate(mip_new = if_else(mip > rad_mip + rad_mip_ci, mip, -mip)) 
pri_data$intertype <- ifelse(pri_data$mro_new > 0 & pri_data$mip_new > 0, "Type 1",ifelse(pri_data$mro_new < 0 & pri_data$mip_new > 0, "Type 2",ifelse(pri_data$mro_new < 0 & pri_data$mip_new < 0, "Type 3","Type 4")))
pri_data$sample <- rownames(pri_data)

del_var <- left_join(del_var,pri_data[,c("intertype","sample")],by = "sample")
del_var <- as.data.frame(del_var)

pop_phyla <- read.table('6555pop.phyla',header = T,row.names = 1)

del_var <- left_join(del_var,pop_phyla,by = 'bin')

del_var_mean <- aggregate(del_var[,3:4],list(del_var$phyla),mean)
phy_pre <- as.data.frame(table(unique(del_var[,c(1,6)])$phyla))
unique(del_var_mean$Group.1 == phy_pre$Var1)
del_var_mean <- data.frame(del_var_mean,phy_pre)
del_var_mean$pre <- del_var_mean$Freq/524

p1 <- ggplot(del_var_mean,aes(mro_var,pre))+
  geom_point(size=3,shape = 17,aes(color = mro_var))+
  scale_color_gradient(low = "#dcdcdc",high = "#3ca6cc",limits = c(-0.001, 0.0002))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=8,colour = 'black'),
        axis.text.y=element_text(size=8,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  labs(x="Mean contributions",y="Prevalence")+
  scale_x_continuous(limits = c(-0.001,0.0002),breaks=seq(-0.001,0.0002,0.0003))+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.25))+
  coord_cartesian(xlim= c(-0.001,0.0002),ylim = c(0,1))+
  expand_limits(x=-0.001,y=0)

p2 <- ggplot(del_var_mean,aes(mip_var,pre))+
  geom_point(size=3,shape = 17,aes(color = mip_var))+
  scale_color_gradient(low = "#dcdcdc",high = "#ff4f42",limits = c(0, 0.6))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=8,colour = 'black'),
        axis.text.y=element_text(size=8,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  labs(x="Mean contributions",y="Prevalence")+
  scale_x_continuous(limits = c(0,0.6),breaks=seq(0,0.6,0.2))+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.25))+
  coord_cartesian(xlim= c(0,0.6),ylim = c(0,1))+
  expand_limits(x=0,y=0)

library(gridExtra)
plots <- list(p1,p2)
grid.arrange(grobs = plots,ncol=1)
