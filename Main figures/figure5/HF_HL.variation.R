library(Hmisc)
library(dplyr)
library(tidyr)
library(lme4)
library(ggsignif)
setwd("D:/Prof.Huang/Global_gut_microbiome/FIGURES4/FIG5")
mro_mip <- read.table('116site.mro_mip_div',sep= "\t",row.names = 1,check.names = F)


sample_data <- read.delim('116sample.data',sep= "\t",row.names = 1,check.names = F)

unique(rownames(mro_mip) == sample_data$`Sample ID`)
sample_data <- sample_data %>%
  mutate(group = case_when(  
    Species == "Nomascus_concolor" & Month %in% c("202104","202107","202201") ~ "HL",  
    Species == "Nomascus_concolor" & !(Month %in% c("202104","202107","202201")) ~ "HF",  
    Species == "Trachypithecus_crepusculus" & Month %in% c("202104","202107","202201") ~ "HL", 
    Species == "Trachypithecus_crepusculus" & !(Month %in% c("202104","202107","202201")) ~ "HF", 
    TRUE ~ NA_character_ # 如果没有匹配任何条件，返回NA  
  ))

total_data <- data.frame(mro_mip,sample_data)
total_data$Sex <- as.factor(total_data$Sex)
total_data$Age <- as.factor(total_data$Age)
total_data$Month <- as.factor(total_data$Month)
total_data$Social.Group <- as.factor(total_data$Social.Group)

total_nc <- total_data[total_data$Species == "Nomascus_concolor",]
total_tc <- total_data[total_data$Species == "Trachypithecus_crepusculus",]

wilcox.test(total_nc$mro~total_nc$group) #higher mro for HL than HF
wilcox.test(total_nc$mip~total_nc$group)
wilcox.test(total_tc$mro~total_tc$group)
wilcox.test(total_tc$mip~total_tc$group)

pairwise.wilcox.test(total_nc$mro,total_nc$Month)
pairwise.wilcox.test(total_nc$mip,total_nc$Month)
pairwise.wilcox.test(total_tc$mro,total_tc$Month)
pairwise.wilcox.test(total_tc$mip,total_tc$Month)

p1 <- ggplot(total_data,aes(group,mro))+
  #geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Diet group",y="Competition (MRO score)")+
  geom_signif(comparisons = list(c("HF","HL")),
              map_signif_level = F, 
              textsize = 3, 
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
  scale_y_continuous(limits = c(0.5,0.65),breaks = seq(0.5,0.65,0.05))+
  coord_cartesian(ylim= c(0.5,0.65))+
  expand_limits(y=0.5)+
  facet_grid(.~Species)

p2 <- ggplot(total_data,aes(group,mip))+
  #geom_jitter(shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Diet group",y="Cooperation (MIP score)")+
  geom_signif(comparisons = list(c("HF","HL")),
              map_signif_level = F, 
              textsize = 3, 
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
  scale_y_continuous(limits = c(100,160),breaks = seq(100,160,20))+
  coord_cartesian(ylim= c(100,160))+
  expand_limits(y=100)+
  facet_grid(.~Species)


#we then plot the community composition of nc and tc between hf and hl
#to identify wherther the varition of diet impact the proportion of baci_a
setwd("D:/Prof.Huang/Global_gut_microbiome/FIGURES4/FIG5")

mat <- read.delim("pop_116sample.abun",header = T,row.names = 1,check.names = F)
phyla <- read.table('6555pop.phyla',header = T)
rownames(phyla) <- phyla$bin
sample_data <- read.delim('116sample.data',header=T, sep="\t", row.names = 1,check.names = F)
rownames(sample_data) <- sample_data$`Sample ID`

unique(rownames(mat) == phyla$bin)
mat <- mat[match(phyla$bin,rownames(mat)),]

phylas <- phyla$phyla
mat_new <- cbind(phylas,mat)

mat_new[,2:ncol(mat_new)] <- lapply(mat_new[,2:ncol(mat_new)],function(x) as.numeric(x != 0))

phy_mat <- aggregate(mat_new[,2:ncol(mat_new)],list(mat_new$phylas),sum)
rownames(phy_mat) <- phy_mat$Group.1
phy_mat <- as.data.frame(t(phy_mat[,-1]))
phy_mat <- phy_mat/rowSums(phy_mat)

unique(rownames(phy_mat) == rownames(sample_data))

stat_phyla <- c("p__Bacillota_A","p__Pseudomonadota","p__Spirochaetota")#,"p__Bacillota","p__Actinomycetota","p__Pseudomonadota")
phy_mat <- phy_mat[,stat_phyla]
#phy_mat$Others <- 1 - rowSums(phy_mat)
sample_data <- sample_data %>%
  mutate(group = case_when(  
    Species == "Nomascus_concolor" & Month %in% c("202104","202107","202201") ~ "HL",  
    Species == "Nomascus_concolor" & !(Month %in% c("202104","202107","202201")) ~ "HF",  
    Species == "Trachypithecus_crepusculus" & Month %in% c("202104","202107","202201") ~ "HL", 
    Species == "Trachypithecus_crepusculus" & !(Month %in% c("202104","202107","202201")) ~ "HF", 
    TRUE ~ NA_character_ # 如果没有匹配任何条件，返回NA  
  ))
spe_diet <- sample_data[,c("Species","group")]
phy_mat <- data.frame(phy_mat,spe_diet)

phy_mat_nc <- phy_mat[phy_mat$Species == "Nomascus_concolor",]
phy_mat_tc <- phy_mat[phy_mat$Species == "Trachypithecus_crepusculus",]
wilcox.test(phy_mat_nc$p__Bacillota_A~phy_mat_nc$group)
wilcox.test(phy_mat_nc$p__Spirochaetota~phy_mat_nc$group)
wilcox.test(phy_mat_nc$p__Pseudomonadota~phy_mat_nc$group)

p3 <- ggplot(phy_mat_nc,aes(group,p__Bacillota_A))+
  #geom_jitter(aes(group,Freq),shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0.1,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Diet group",y="Proportion of Bacillota_A")+
  geom_signif(comparisons = list(c("HF","HL")),
              map_signif_level = F, 
              textsize = 3, 
              step_increase = 0.1,
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.3,0.6),breaks = seq(0.3,0.6,0.1))+
  coord_cartesian(ylim= c(0.3,0.6))+
  expand_limits(y=0.3)

ps2 <- ggplot(phy_mat_tc,aes(group,p__Bacillota_A))+
  #geom_jitter(aes(group,Freq),shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0.1,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Diet group",y="Proportion of Bacillota_A")+
  geom_signif(comparisons = list(c("HF","HL")),
              map_signif_level = F, 
              textsize = 3, 
              step_increase = 0.1,
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.7,0.9),breaks = seq(0.7,0.9,0.1))+
  coord_cartesian(ylim= c(0.7,0.9))+
  expand_limits(y=0.7)

#we then want to know whether the increase in the proportion of baci_a was resulted from the baci_a with GH5
setwd("D:/Prof.Huang/Global_gut_microbiome/FIGURES4/Fig.S8")
cazy <- read.delim('6555genomes.pro.ko_ano_addphyla_cazy.baci_a.cut.rmr',header =F)
colnames(cazy) <- c("phyla","bin","pc","ko","fun","fam")
baci_a_gh5 <- unique(cazy[cazy$phyla == "p__Bacillota_A" & cazy$fam == "GH5",]$bin)

baci_a_gh5_div <- colSums(mat_new[baci_a_gh5,2:ncol(mat_new)])
baci_a_gh5_prop <- baci_a_gh5_div/colSums(mat_new[,2:ncol(mat_new)])
phy_mat2 <- data.frame(baci_a_gh5_prop,spe_diet)
phy_mat2_nc <- phy_mat2[phy_mat2$Species == "Nomascus_concolor",]
wilcox.test(phy_mat2_nc$baci_a_gh5_prop~phy_mat_nc$group)

p4 <- ggplot(phy_mat2_nc,aes(group,baci_a_gh5_prop))+
  #geom_jitter(aes(group,Freq),shape = 16,size=1.5,width = 0.15,alpha = 0.5,color = "#c8c8c8")+
  stat_boxplot(geom = "errorbar",width=0.1,aes(color = group),position = position_dodge(1),alpha = 0)+
  geom_boxplot(aes(color = group),width=0.3,position = position_dodge(1),outlier.shape = NA,fill = "#FFFFFF",alpha = 0)+
  scale_color_manual(values=c("#ff4f42","#3ca6cc"))+
  labs(x = "Diet group",y="Proportion of Bacillota_A with GH5")+
  geom_signif(comparisons = list(c("HF","HL")),
              map_signif_level = F, 
              textsize = 3, 
              step_increase = 0.1,
              test = "wilcox.test")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=8),
        legend.position = 'top')+
  #scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0.2,0.5),breaks = seq(0.2,0.5,0.1))+
  coord_cartesian(ylim= c(0.2,0.5))+
  expand_limits(y=0.2)

library(gridExtra)
plots <- list(p1,p2,p3,p4)
grid.arrange(grobs = plots,ncol=4)