pop_char <- read.table('6555pop.char.table',header = T)

pop_char_good <- pop_char[pop_char$quality > 95,]

p1 <- ggplot(pop_char_good,aes(in_num,out_num))+
  geom_point(size=3,shape = 16,aes(color = mro))+
  scale_color_gradient(low = "#f0f0f0",high = "#3ca6cc",limits = c(0.3, 0.7))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  labs(x="Number of input compounds",y="Number of output compounds")+
  scale_x_continuous(limits = c(40,280),breaks = seq(40,280,80))+
  scale_y_continuous(limits = c(0,60),breaks=seq(0,60,20))+
  coord_cartesian(xlim= c(40,280),ylim = c(0,60))+
  expand_limits(x=40,y=0) 

p2 <- ggplot(pop_char_good,aes(in_num,out_num))+
  geom_point(size=3,shape = 16,aes(color = mip))+
  scale_color_gradient(low = "#f0f0f0",high = "#ff4f42",limits = c(0, 5))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black',angle = 90),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  labs(x="Number of input compounds",y="Number of output compounds")+
  scale_x_continuous(limits = c(40,280),breaks = seq(40,280,80))+
  scale_y_continuous(limits = c(0,60),breaks=seq(0,60,20))+
  coord_cartesian(xlim= c(40,280),ylim = c(0,60))+
  expand_limits(x=40,y=0) 

pop_char_good$phyla2 <- ifelse(pop_char_good$phyla %in% c("p__Bacillota_A","p__Pseudomonadota"),pop_char_good$phyla,"Others")
pop_char_good$phyla2 <- factor(pop_char_good$phyla2,levels = c("p__Bacillota_A","p__Pseudomonadota","Others"))
pop_char_good$size <- pop_char_good$length/1000
p3 <- ggplot(pop_char_good)+
  geom_point(aes(size,in_num,color = phyla2),size=3)+
  geom_smooth(aes(size,in_num),method='lm',fill=NA,color = "#000000")+
  scale_color_manual(values=c("#3ca6cc","#ff4f42","#bebebe"))+
  labs(x="Genome size (kb)",y="Number of input compounds")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  expand_limits(x=0,y=40)+
  scale_x_continuous(limits = c(0,9000),breaks = seq(0,9000,3000))+
  scale_y_continuous(limits = c(40,280),breaks = seq(40,280,80))

p4 <- ggplot(pop_char_good)+
  geom_point(aes(size,out_num,color = phyla2),size=3)+
  geom_smooth(aes(size,out_num),method='lm',fill=NA,color = "#000000")+
  scale_color_manual(values=c("#3ca6cc","#ff4f42","#bebebe"))+
  labs(x="Genome size (kb)",y="Number of output compounds")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'none')+
  expand_limits(x=0,y=0)+
  scale_x_continuous(limits = c(0,9000),breaks = seq(0,9000,3000))+
  scale_y_continuous(limits = c(0,60),breaks=seq(0,60,20))

library(gridExtra)
plots <- list(p1,p2,p3,p4)
grid.arrange(grobs = plots,ncol=2)
