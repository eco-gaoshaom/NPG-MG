all_curve <- read.table("npgmg.curve_data",header = T)
nons_curve <- read.table("npgmg.nonsingleton.curve_data",header = T)
all_curve <- rbind(c(0,0,0),all_curve)
nons_curve <- rbind(c(0,0,0),nons_curve)

curve_data <- rbind(all_curve,nons_curve)
curve_data$group <- c(replicate(nrow(all_curve),"All species"),replicate(nrow(nons_curve),"No singletons"))
curve_data$group <- factor(curve_data$group,levels = c("All species","No singletons"))

ggplot(curve_data)+
  geom_line(aes(x=xaxis,y=X1,color= group),linewidth=2)+
  scale_color_manual(values=c(rgb(104,102,135,max=255),rgb(220,220,220,max=255)))+
  theme_bw()+
  labs(x = "Number of genomes", y = "Number of species")+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth=0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'right')+expand_limits(x=0,y=0)+
  scale_x_continuous(limits =c(0,40000),breaks = seq(0,40000,8000))+
  scale_y_continuous(limits =c(0,7500),breaks = seq(0,7500,1500))+
  coord_cartesian(ylim=c(0,7500),xlim = c(0,40000))+
  geom_hline(yintercept=max(all_curve$X1),linetype=4,color="grey",linewidth=1)+ 
  geom_vline(xintercept=max(all_curve$xaxis),linetype=4,color="grey",linewidth=1)+
  geom_hline(yintercept=max(nons_curve$X1),linetype=4,color="grey",linewidth=1)+ 
  geom_vline(xintercept=max(nons_curve$xaxis),linetype=4,color="grey",linewidth=1)
