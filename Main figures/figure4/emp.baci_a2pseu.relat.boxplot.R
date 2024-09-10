library(ggplot2)
library(ggsignif)
phy <- read.table('emp_deblur_150bp.release1.table.sed.gtdb_phy_level',header = T,check.names = F)
colnames(phy) <- c("phy",colnames(phy)[1:(ncol(phy)-1)])
rownames(phy) <- phy[,1]
phy <- phy[,-1]
phy_rela <- t(phy)/as.numeric(rowSums(t(phy)))
rownames(phy_rela) <- colnames(phy)

meta <- read.delim('emp_qiime_mapping_release1_20170912.tsv',header = T)
rownames(meta) <- meta$X.SampleID
phy_meta <- meta[rownames(phy_rela),]

phy_rela <- as.data.frame(phy_rela)
phy_rela$group <- phy_meta$empo_1
group_rela <- aggregate(phy_rela[,c("Bacillota_A","Pseudomonadota")],list(phy_rela$group),mean)

phy_rela <- phy_rela[phy_rela$group != "Control",c("Bacillota_A","Pseudomonadota","group")]
phy_rela_df <- reshape2::melt(phy_rela,id.vars = c("group"))

wilcox.test(phy_rela$Pseudomonadota~phy_rela$group)
wilcox.test(phy_rela$Bacillota_A~phy_rela$group)

phy_rela_df$group <- factor(phy_rela_df$group,levels = c("Host-associated","Free-living"))

ggplot(phy_rela_df,aes(variable,value,fill=group))+
  #geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Phyla",y="Relative abundance")+
  geom_signif(comparisons = list(c("Host-associated","Free-living")),
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
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))+
  coord_cartesian(ylim= c(0,1))+
  expand_limits(y=0)
