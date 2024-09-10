#we would like to use the three phyla to predict intertypes
library(dplyr)
library(randomForest)
library(caret)
library(vegan)

pri_data <- read.delim("site.rand_mro_mip_div",header = T,check.names = F,row.names = 1)
pri_data <- pri_data %>%  mutate(mro_new = if_else(mro > rad_mro + rad_mro_ci, mro, -mro)) 
pri_data <- pri_data %>%  mutate(mip_new = if_else(mip > rad_mip + rad_mip_ci, mip, -mip)) 
pri_data$intertype <- ifelse(pri_data$mro_new > 0 & pri_data$mip_new > 0, "Type 1",ifelse(pri_data$mro_new < 0 & pri_data$mip_new > 0, "Type 2",ifelse(pri_data$mro_new < 0 & pri_data$mip_new < 0, "Type 3","Type 4")))

pri_data$sample <- rownames(pri_data)
pri_db <- pri_data[,c("sample","intertype")]

popu_phy <- read.table("6555pop.phyla",header = T)
phy_div <- as.data.frame(table(popu_phy$phyla))
phy_div <- phy_div[order(phy_div$Freq,decreasing = T),]
good_phy <- as.character(phy_div[1:6,]$Var1)

data <- read.table("npgmg_524subsampled.phyla_abun",header = T,check.names = F)
data[,2:ncol(data)] <- lapply(data[,2:ncol(data)],function(x) as.numeric(x != 0))

new <- aggregate(data[,2:ncol(data)],list(data$phyla),sum)
rownames(new) <- new$Group.1
new <- as.data.frame(t(new[,-1]))

new$sample <- rownames(new)

unique(rownames(new) == rownames(pri_db))

new <- dplyr::left_join(new,pri_db,by = "sample")
new <- new[, -which(names(new) == "sample")]
rownames(new) <- rownames(pri_data)

new1 <- decostand(new[,good_phy],method = "standardize")

df <- data.frame(new1,new$intertype)
colnames(df) <- c(colnames(new1),"intertype")

#since the number of the four intertypes are unbalanced, we use smote to create new data
library(smotefamily)

imp_list <- list()
mda_mdg <- list()
accu_list <- list()

sample_with_all_types <- function(df, total_samples) {  
  while (TRUE) {  
    set.seed(sample.int(.Machine$integer.max, 1)) 
    sampled_df <- df %>%  
      sample_n(size = total_samples)
    if (length(unique(sampled_df$intertype)) == length(unique(df$intertype)) &
        min(table(sampled_df$intertype)) >= 5) {  
      return(sampled_df)  
    }  
  }  
}

for (i in 1:999){
  print(i)
  df_tests <- sample_with_all_types(df,131)
  df_trains <- df[setdiff(rownames(df),rownames(df_tests)),]
  
  df_trains_t1 <- df_trains[df_trains$intertype == "Type 1",]
  df_trains_t2 <- df_trains[df_trains$intertype == "Type 2",]
  df_trains_t3 <- df_trains[df_trains$intertype == "Type 3",]
  df_trains_t4 <- df_trains[df_trains$intertype == "Type 4",]
  
  df_trains_t4_sm <- SMOTE(df_trains_t4[,-ncol(df_trains_t1)],df_trains_t4[,ncol(df_trains_t1)],K = 5,dup_size = 3)$data
  df_trains_t1_sm <- SMOTE(df_trains_t1[,-ncol(df_trains_t1)],df_trains_t1[,ncol(df_trains_t1)],K = 5,dup_size = floor(nrow(df_trains_t4_sm)/nrow(df_trains_t1) - 1))$data
  df_trains_t2_sm <- SMOTE(df_trains_t2[,-ncol(df_trains_t1)],df_trains_t2[,ncol(df_trains_t1)],K = 5,dup_size = floor(nrow(df_trains_t4_sm)/nrow(df_trains_t2) - 1))$data
  df_trains_t3_sm <- SMOTE(df_trains_t3[,-ncol(df_trains_t1)],df_trains_t3[,ncol(df_trains_t1)],K = 5,dup_size = floor(nrow(df_trains_t4_sm)/nrow(df_trains_t3) - 1))$data

  #colnames(df_trains_t4) <- colnames(df_trains_t1_sm)
  df_trains_df <- rbind(df_trains_t1_sm,df_trains_t2_sm,df_trains_t3_sm,df_trains_t4_sm)
  
  df_trains_df$class <- as.factor(df_trains_df$class)
  rf_balanced <- randomForest(class ~ ., data = df_trains_df, ntree = 1000, mtry = 3, importance = TRUE)  
  imp <- importance(rf_balanced)
  imp_list[[i]] <- imp[,1:4]
  mda_mdg[[i]] <- imp[,5:6]
  accuracy0 <- sum(rf_balanced$predicted == df_trains_df$class)/nrow(df_trains_df)
  
  train_acc <- c(accuracy0,1-rf_balanced$confusion[,5])
  
  df_tests_pd <- predict(rf_balanced,df_tests)
  df_tests$intertype <- as.factor(df_tests$intertype)
  confu_mats <- confusionMatrix(df_tests_pd,df_tests$intertype)
  test_acc <- c(confu_mats$overall[1],confu_mats$byClass[,11])
  accu_list[[i]] <- data.frame(train_acc,test_acc)
}
setwd("D:/Prof.Huang/Global_gut_microbiome/FIGURES4/FIG3")
save.image("random_forest.RData")

load("random_forest.RData")
imp_mean <- as.data.frame(sapply(1:4, function(j) Reduce(`+`, lapply(imp_list, function(df) df[, j]))) / length(imp_list))
colnames(imp_mean) <- c("Type 1","Type 2","Type 3","Type 4")

mda_mdg_mean <- as.data.frame(sapply(1:2, function(j) Reduce(`+`, lapply(mda_mdg, function(df) df[, j]))) / length(mda_mdg))
colnames(mda_mdg_mean) <- c("MDA","MDG")

accu_mean <- as.data.frame(sapply(1:2, function(j) Reduce(`+`, lapply(accu_list, function(df) df[, j]))) / length(accu_list))
colnames(accu_mean) <- c("Train","Test")
rownames(accu_mean) <- c("Total","Type 1","Type 2","Type 3","Type 4")

#plot heatmap
imp_mean$phyla <- rownames(imp_mean)
imp_mean_melt <- reshape2::melt(imp_mean,id.vars = c("phyla"))
imp_mean_melt$phyla <- factor(imp_mean_melt$phyla,levels = rownames(imp_mean))
imp_mean_melt$variable <- factor(imp_mean_melt$variable,levels = colnames(imp_mean)[1:4])

ggplot(imp_mean_melt,aes(variable,phyla))+
  geom_point(aes(color =value),size = 8,alpha=1,shape=15)+
  scale_size_area(max_size = 12,breaks = 2)+
  theme_bw()+coord_flip()+
  scale_color_gradient(low = "#ffffff",high = "#ab84b6",limits = c(40, 160))+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=8,angle= 45, vjust= 1,hjust = 1),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

#plot mda and mdg barplot
mda_mdg_mean$phyla <- rownames(mda_mdg_mean)
mda_mdg_mean_melt <- reshape2::melt(mda_mdg_mean,id.vars = c("phyla"))
mda_mdg_mean_melt$phyla <- factor(mda_mdg_mean_melt$phyla,levels = rownames(mda_mdg_mean))
mda_mdg_mean_melt$variable <- factor(mda_mdg_mean_melt$variable,levels = colnames(mda_mdg_mean)[1:2])
ggplot(mda_mdg_mean_melt,aes(phyla,value,fill=variable))+
  geom_bar(stat = 'identity',position = 'dodge',width =0.3)+
  scale_fill_manual(values = c("#905797","#9b87be"))+
  labs(x="Phyla",y='Value')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(size=1,colour="black"),legend.text = element_text(size=13),
        legend.position = 'right')+
  scale_y_continuous(breaks = seq(0,800,400))+
  coord_cartesian(ylim = c(0,800))

#accuracy barplot
accu_mean$types <- rownames(accu_mean)
accu_mean_melt <- reshape2::melt(accu_mean,id.vars = c("types"))
accu_mean_melt$types <- factor(accu_mean_melt$types,levels = rownames(accu_mean))
accu_mean_melt$variable <- factor(accu_mean_melt$variable,levels = colnames(accu_mean)[1:2])
ggplot(accu_mean_melt,aes(types,value,fill=variable))+
  geom_bar(stat = 'identity',position = 'dodge',width =0.3)+
  scale_fill_manual(values = c("#5bb288","#90bb97"))+
  labs(x="Phyla",y='Accuracy')+
  theme_bw()+coord_flip(ylim = c(0,1))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(size=1,colour="black"),legend.text = element_text(size=13),
        legend.position = 'right')+
  scale_y_continuous(breaks = seq(0,1,0.5))
  
#plot good_phy proportion in four intertypes
new_del <- new[,-ncol(new)]
new_prop <- new_del/rowSums(new_del)
prop <- new_prop[,good_phy]
prop$Others <- 1 - rowSums(prop)
prop$intertype <- new$intertype
prop_mean <- aggregate(prop[,1:7],list(prop$intertype),mean)

prop_melt <- reshape2::melt(prop_mean,id.vars = c("Group.1"))
prop_melt$Group.1 <- factor(prop_melt$Group.1,levels = c("Type 1","Type 2","Type 3","Type 4"))
prop_melt$variable <- factor(prop_melt$variable,levels = c(good_phy,"Others"))

ggplot(prop_melt)+
  geom_bar(aes(Group.1,value,fill=variable),
           stat = 'identity',position = 'stack',width =0.4)+
  scale_fill_manual(values = c("#3ca6cc","#ff4f42","#5bb288","#ab84b6","#FAC00F","#BC587D","#dcdcdc"))+
  theme_bw()+coord_flip()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=5,colour = 'black',angle = 90,hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=6),
        legend.position = 'top')+
  labs(x="Intertypes",y="Proportion of populations")+
  #scale_x_continuous(limits = c(0.45,0.6),breaks=seq(0.45,0.6,0.05))+
  scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.25))+
  expand_limits(y=0)
