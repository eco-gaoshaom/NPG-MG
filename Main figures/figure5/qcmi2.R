library("qcmi")
library("phyloseq")
library("ggClusterNet")
library("tidyverse")
library("WGCNA")
library("igraph")
library("Matrix")
library("vegan")
library("magrittr")
library("reticulate")
library("SpiecEasi")
library("ggstatsplot")
library("dplyr")

mat <- read.delim("pop_116sample.abun",header = T,row.names = 1,check.names = F)

sample_data <- read.delim('116sample.data',header=T, sep="\t", row.names = 1,check.names = F)
rownames(sample_data) <- sample_data$`Sample ID`
sample_data <- sample_data %>%
  mutate(group = case_when(  
    Species == "Nomascus_concolor" & Month %in% c("202104","202107","202201") ~ "HL",  
    Species == "Nomascus_concolor" & !(Month %in% c("202104","202107","202201")) ~ "HF",  
    Species == "Trachypithecus_crepusculus" & Month %in% c("202104","202107","202201") ~ "HL", 
    Species == "Trachypithecus_crepusculus" & !(Month %in% c("202104","202107","202201")) ~ "HF", 
    TRUE ~ NA_character_ # 如果没有匹配任何条件，返回NA  
  ))

unique(colnames(mat) == rownames(sample_data))

mat_nc <- mat[,1:68]
mat_tc <- mat[,69:116]
mat_nc <- mat_nc[rowSums(mat_nc) > 0,]
mat_tc <- mat_tc[rowSums(mat_tc) > 0,]

phyla <- read.table('6555pop.phyla',header = T)
rownames(phyla) <- phyla$bin
phyla <- phyla[,c(1,2)]
tax_nc <- phyla[rownames(mat_nc),]
tax_tc <- phyla[rownames(mat_tc),]

date_nc <- as.data.frame(as.Date(substr(colnames(mat_nc), 1, 6), format = "%y%m%d"))
date_tc <- as.data.frame(as.Date(substr(colnames(mat_tc), 2, 7), format = "%y%m%d"))
rownames(date_nc) <- colnames(mat_nc)
rownames(date_tc) <- colnames(mat_tc)
colnames(date_nc) <- colnames(date_tc) <- c("date")

save.wd="./result2"
if(!dir.exists(save.wd)){dir.create(save.wd)}
setwd(save.wd)

ps_nc <- trans_ps(otu_table = mat_nc, taxa_table =  tax_nc)
print(ps_nc)

ps_tc <- trans_ps(otu_table = mat_tc, taxa_table =  tax_tc)
print(ps_tc)

filteredps_nc <- filter_ps(ps = ps_nc, occurrence.threshold=17,  abundance.threshold=502)
filteredps_tc <- filter_ps(ps = ps_tc, occurrence.threshold=12,  abundance.threshold=495)
print(filteredps_nc)
print(filteredps_tc)

filteredps_nc = filter_OTU_ps(ps = filteredps_nc,Top = 100)
filteredps_tc = filter_OTU_ps(ps = filteredps_tc,Top = 257)
save.image('qcmi.1.RData')


library("phyloseq")
library("SpiecEasi")
load("qcmi.1.RData")
cal_net <- function(ps = filteredps, N = 100,r.threshold=0.6,p.threshold=0.05,
                    lambda.min.ratio=1e-2,nlambda=20,ncores=10,
                    method ="spieceasi"){
  method = "spieceasi"
  se.net <- spiec.easi(ps, method='mb', lambda.min.ratio=0.01, nlambda=20, pulsar.params=list(rep.num=50, ncores=10))
  ig.mb  <- adj2igraph(getRefit(se.net), vertex.attr= list(label = row.names(ps@otu_table)))
  
  sebeta <- symBeta(getOptBeta(se.net), mode='maxabs')
  elist.mb <- summary(sebeta)
  
  
  result=list(sebeta,elist.mb,method,ps,ig.mb,se.net)
  names(result)[1] <- "occor.r"
  names(result)[2] <- "elist"
  names(result)[3] <- "method"
  names(result)[4] <- "phyloseq"
  names(result)[5] <- "igraph"
  names(result)[6] <- "result.se"
  return(result)
}

result_net_nc <- cal_net(ps = filteredps_nc,N=100,ncores=10,lambda.min.ratio=1e-2,nlambda=20,method ="spieceasi")
summary(result_net_nc)

result_net_tc <-  cal_net(ps = filteredps_tc,N=257,ncores=10,lambda.min.ratio=1e-2,nlambda=20,method ="spieceasi")
summary(result_net_tc)
#write.csv(result_net_nc$elist,file="nc.edgelist.csv")
save.image("qcmi.2.RData")


load("qcmi.2.RData")
result_net_nc$elist <- summary(result_net_nc$occor.r)
result_net_tc$elist <- summary(result_net_tc$occor.r)

ps_sub_nc = filter_OTU_ps(ps = filteredps_nc,Top = 100)
otu_table_nc_raw = as.data.frame(t(vegan_otu(ps_sub_nc)))
otu_table_nc = as.data.frame(t(vegan_otu(ps_sub_nc)))
row.names(otu_table_nc)<-c(1:100)

ps_sub_tc = filter_OTU_ps(ps = filteredps_tc,Top = 257)
otu_table_tc_raw = as.data.frame(t(vegan_otu(ps_sub_tc)))
otu_table_tc = as.data.frame(t(vegan_otu(ps_sub_tc)))
row.names(otu_table_tc)<-c(1:257)


TEST_link_dl = function(link_table_row, OTUabd, geodist){
  OTU1 = as.character(as.matrix(link_table_row[1]))
  OTU2	 = as.character(as.matrix(link_table_row[2]))
  OTU1abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU1),])
  OTU2abd = as.numeric(OTUabd[which(row.names(OTUabd)== OTU2),])
  distOTU1 = vegdist(cbind(OTU1abd,rep(0,length(OTU1abd))), method="bray")
  distOTU2 = vegdist(cbind(OTU2abd,rep(0,length(OTU2abd))), method="bray")
  
  cor1 = cor.test(distOTU1, geodist)
  cor2 = cor.test(distOTU2, geodist)
  
  
  return(c(OTU1, OTU2, cor1$estimate, cor1$p.value, cor2$estimate, cor2$p.value))
}

assigned_proc = function(link_table_row, OTUabd, p=0.05 ,data,cutoff=0,method=c("dl", "ef")){
  
  if (method %in% c("dl")) {
    
    edges =link_table_row
    OTUs = OTUabd 
    distances = data
    
    OTUs_order = as.data.frame(OTUs[,order(names(OTUs))])
    distances_order = as.data.frame(distances[order(row.names(distances)),])
    row.names(distances_order) <- colnames(OTUs_order)
    geodist = dist(distances_order[,1])
    
    test_result = t(apply(edges, 1, FUN=test_link_dl, OTUabd=OTUs_order, geodist=geodist))
    test_result=as.data.frame(test_result)
    #test_result=cbind(row.names(test_result),test_result)
    colnames(test_result)=c("OTU1","OTU2","cor1","P1","cor2","P2")
    #scene1: positive links   double positive covary with geo
    test_result_pick1 = test_result[which((test_result['cor1']>cutoff) & (test_result['cor2']>cutoff) & (test_result['P1']<=p) & (test_result['P2']<=p)),]
    
    #scene2: positive links   double negative covary with geo
    test_result_pick2 = test_result[which((test_result[,3]<cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
    
    
    result= rbind(test_result_pick1 ,test_result_pick2)
    result_dl= result
    
    dl=rep("yes",times=nrow(result_dl))
    result_dl=cbind(result_dl,dl)
    
    return(result_dl)
  } 
  
  if (method %in% c("ef")) {
    
    edges =link_table_row
    OTUs = OTUabd 
    distances = data
    
    OTUs_order = as.data.frame(OTUs[,order(names(OTUs))])
    distances_order = as.data.frame(distances[order(row.names(distances)),])
    envdist = vegdist(distances_order, method="bray")
    
    
    #scene1: positive links   double positive covary with env
    edges1=edges[which((edges[,3]>=cutoff) ),]
    test_result = t(apply(edges1, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
    test_result=as.data.frame(test_result)
    test_result_pick1 = test_result[which((test_result[,3]>cutoff) & (test_result[,5]>cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
    
    #scene2: positive links   double negative covary with env
    edges2=edges[which((edges[,3]>=cutoff) ),]
    test_result = t(apply(edges2, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
    test_result=as.data.frame(test_result)
    test_result_pick2 = test_result[which((test_result[,3]<cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
    
    #scene3: negative links   1- 1+  covary with env
    edges3=edges[which((edges[,3]<cutoff) ),]
    test_result = t(apply(edges3, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
    test_result=as.data.frame(test_result)
    test_result_pick3 = test_result[which((test_result[,3]<cutoff) & (test_result[,5] > cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
    
    #scene4: negative links   1+ 1- covary with env
    edges4=edges[which((edges[,3]<cutoff) ),]
    test_result = t(apply(edges4, 1, FUN=test_link_env, OTUabd=OTUs_order, envdist=envdist))
    test_result=as.data.frame(test_result)
    test_result_pick4 = test_result[which((test_result[,3]>cutoff) & (test_result[,5]<cutoff) & (test_result[,4]<=p) & (test_result[,6]<=p)),]
    
    
    result= rbind(test_result_pick1 ,test_result_pick2,test_result_pick3,test_result_pick4)
    
    
    result_ef= result
    
    colnames(result_ef)=c("OTU1","OTU2","cor1","P1","cor2","P2")
    
    ef=rep("yes",times=nrow(result_ef))
    
    result_ef=cbind(result_ef,ef)
    
    return(result_ef)
  } 
}

result_nc_dl <- assigned_proc(link_table_row = result_net_nc$elist,OTUabd = otu_table_nc,p=0.05 , data= date_nc,cutoff=0, method="dl")
result_tc_dl <- assigned_proc(link_table_row = result_net_tc$elist,OTUabd = otu_table_tc,p=0.05 , data= date_tc,cutoff=0, method="dl")


library(Hmisc)
diet_nc <- sample_data[1:68,9:12]
diet_tc <- sample_data[69:116,9:12]

plot(varclus(as.matrix(diet_nc)))
result_nc_leaf <- assigned_process(link_table_row=result_net_nc$elist, OTUabd=otu_table_nc, p=0.05 , data= diet_nc['Leaf'],cutoff=0, method="ef")
result_tc_leaf <- assigned_process(link_table_row=result_net_tc$elist, OTUabd=otu_table_tc, p=0.05 , data= diet_tc['Leaf'],cutoff=0, method="ef")
result_nc_flower <- assigned_process(link_table_row=result_net_nc$elist, OTUabd=otu_table_nc, p=0.05 , data= diet_nc['Flower'],cutoff=0, method="ef")
result_tc_flower <- assigned_process(link_table_row=result_net_tc$elist, OTUabd=otu_table_tc, p=0.05 , data= diet_tc['Flower'],cutoff=0, method="ef")

library(dplyr)
total_link_nc=row.names(result_net_nc$elist)
dl_link_nc=as.character(row.names(result_nc_dl))
ef_link_nc=union(as.character(row.names(result_nc_leaf)),as.character(row.names(result_nc_flower)))
bi_link_nc = setdiff(total_link_nc,union(dl_link_nc,ef_link_nc))
biedge_nc = as.data.frame(result_net_nc$elist[bi_link_nc,])

total_link_tc=row.names(result_net_tc$elist)
dl_link_tc=as.character(row.names(result_tc_dl))
ef_link_tc=union(as.character(row.names(result_tc_leaf)),as.character(row.names(result_tc_flower)))
bi_link_tc = setdiff(total_link_tc,union(dl_link_tc,ef_link_tc))
biedge_tc = as.data.frame(result_net_tc$elist[bi_link_tc,])


nc_100tax <- phyla[rownames(otu_table_nc_raw),]
tc_257tax <- phyla[rownames(otu_table_tc_raw),]

nc_100tax$bin <- c(1:100)
colnames(nc_100tax) <- c("i","phyla")
biedge_nc <- left_join(biedge_nc,nc_100tax,by="i")
colnames(nc_100tax) <- c("j","phyla")
biedge_nc <- left_join(biedge_nc,nc_100tax,by="j")

nc_100tax$phyla <- rownames(nc_100tax)
colnames(nc_100tax) <- c("i","bin")
biedge_nc <- left_join(biedge_nc,nc_100tax,by = "i")
colnames(nc_100tax) <- c("j","bin")
biedge_nc <- left_join(biedge_nc,nc_100tax,by = "j")

tc_257tax$bin <- c(1:257)
colnames(tc_257tax) <- c("i","phyla")
biedge_tc <- left_join(biedge_tc,tc_257tax,by="i")
colnames(tc_257tax) <- c("j","phyla")
biedge_tc <- left_join(biedge_tc,tc_257tax,by="j")

tc_257tax$phyla <- rownames(tc_257tax)
colnames(tc_257tax) <- c("i","bin")
biedge_tc <- left_join(biedge_tc,tc_257tax,by = "i")
colnames(tc_257tax) <- c("j","bin")
biedge_tc <- left_join(biedge_tc,tc_257tax,by = "j")

write.table(biedge_nc,"nc.net.bio_edge",sep = "\t")
write.table(biedge_tc,"tc.net.bio_edge",sep = "\t")

biedge_nc <- read.table("nc.net.bio_edge",header = T)
biedge_tc <- read.table("tc.net.bio_edge",header = T)

#we firstly assessed the mro/mip for adj net, to assess whether the mro/mip correlated with biotic correlations
setwd("D:/Prof.Huang/Global_gut_microbiome/FIGURES3/FIG3")
library(data.table)
xmlname <- read.table("6555_sub.xmlname",header = F)
mro_mat <- fread("6555.mro.index.matrix",header = F)
mip_mat <- fread("6555.mip.index.matrix",header = F)

mro_mat <- as.matrix(mro_mat[,-1])
mip_mat <- as.matrix(mip_mat[,-1])

colnames(mro_mat) <- rownames(mro_mat) <- xmlname$V1
colnames(mip_mat) <- rownames(mip_mat) <- xmlname$V1

nc_pairs <- c()
mro_nc <- c()
mip_nc <- c()
for (i in 1:nrow(biedge_nc)){
  nc_pairs[i*2 - 1] <- paste(biedge_nc[i,6],biedge_nc[i,7])
  nc_pairs[i*2] <- paste(biedge_nc[i,7],biedge_nc[i,6])
  mro_nc[i] <- (mro_mat[biedge_nc[i,6],biedge_nc[i,7]] + mro_mat[biedge_nc[i,7],biedge_nc[i,6]])
  mip_nc[i] <- (mip_mat[biedge_nc[i,6],biedge_nc[i,7]] + mip_mat[biedge_nc[i,7],biedge_nc[i,6]])
}

tc_pairs <- c()
mro_tc <- c()
mip_tc <- c()
for (i in 1:nrow(biedge_tc)){
  tc_pairs[i*2 - 1] <- paste(biedge_tc[i,6],biedge_tc[i,7])
  tc_pairs[i*2] <- paste(biedge_tc[i,7],biedge_tc[i,6])
  mro_tc[i] <- (mro_mat[biedge_tc[i,6],biedge_tc[i,7]] + mro_mat[biedge_tc[i,7],biedge_tc[i,6]])
  mip_tc[i] <- (mip_mat[biedge_tc[i,6],biedge_tc[i,7]] + mip_mat[biedge_tc[i,7],biedge_tc[i,6]])
}

cor.test(biedge_nc$x,mro_nc)
cor.test(biedge_nc$x,mip_nc)

cor.test(biedge_tc$x,mro_tc)
cor.test(biedge_tc$x,mip_tc)

nc_df <- data.frame(biedge_nc,mro_nc,mip_nc)
tc_df <- data.frame(biedge_tc,mro_tc,mip_tc)

write.table(nc_df,"nc.net_edge.mro_mip",sep = "\t")
write.table(tc_df,"tc.net_edge.mro_mip",sep = "\t")

nc_df <- read.table("nc.net.bio_edge",header = T)
tc_df <- read.table("tc.net.bio_edge",header = T)

p5 <- ggplot(nc_df)+
  geom_point(aes(mro_nc,x),size=2.5,color="#969696")+
  geom_smooth(aes(mro_nc,x),method='lm',fill=NA,color = "#3ca6cc")+
  labs(x="MRO score",y="Correlations")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=0.3,y=-0.3)+
  scale_x_continuous(limits = c(0.3,0.9),breaks = seq(0.3,0.9,0.3))+
  scale_y_continuous(limits = c(-0.3,0.3),breaks = seq(-0.3,0.3,0.3))



p6 <- ggplot(tc_df)+
  geom_point(aes(mro_tc,x),size=2.5,color="#969696")+
  geom_smooth(aes(mro_tc,x),method='lm',fill=NA,color = "#3ca6cc")+
  labs(x="MRO score",y="Correlations")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,angle = 90,colour = 'black',hjust = 0.5),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'bottom')+
  expand_limits(x=0.2,y=-0.4)+
  scale_x_continuous(limits = c(0.2,1),breaks = seq(0.2,1,0.4))+
  scale_y_continuous(limits = c(-0.4,0.6),breaks=seq(-0.4,0.6,0.5))

library(gridExtra)
plots <- list(p5,p6)
grid.arrange(grobs = plots,ncol=2)