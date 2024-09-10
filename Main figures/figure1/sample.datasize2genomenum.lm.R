data <- read.delim('803sample_info.txt',header = T)
sample_size <- data[,c("fq_name","reads_num")]

data <- read.delim('39622genome_info.txt',header = T)
sample_genome <- unique(data[,c("genome","sample")])
sample_genome_num <- aggregate(sample_genome$genome,list(sample_genome$sample),length)
sample_pop <- unique(data[,c("sample","rep_pop")])
sample_pop_num <- aggregate(sample_pop$rep_pop,list(sample_pop$sample),length)


overlap <- intersect(sample_size$fq_name,sample_pop_num$Group.1)

sample_size_overlap <- sample_size[sample_size$fq_name %in% overlap,]
sample_pop_num_overlap <- sample_pop_num[sample_pop_num$Group.1 %in% overlap,]

sample_size_overlap <- sample_size_overlap[order(sample_size_overlap$fq_name,decreasing = F),]
sample_pop_num_overlap <- sample_pop_num_overlap[order(sample_pop_num_overlap$Group.1,decreasing = F),]


size2num <- data.frame(sample_size_overlap,sample_pop_num_overlap)
size2num$scale <- (size2num$reads_num)/10^6

summary(lm(size2num$x~size2num$scale))
library(ggplot2)
ggplot(size2num)+
  geom_point(aes(scale,x),size=2,shape=16)+
  geom_smooth(aes(scale,x),method='lm',fill=NA)+
  labs(x="Number of reads (x 10^6)",y="Number of genomes")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust = 0.5,angle = 0,size=8,colour = 'black'),
        axis.text.y=element_text(size=15,colour = 'black'),panel.border = element_blank(),
        axis.line = element_line(linewidth =0.5,colour="black"),legend.text = element_text(size=13),
        legend.position = 'right')+
  scale_x_continuous(limits = c(0,400),breaks = seq(0,400,100))+
  scale_y_continuous(limits = c(0,250),breaks = seq(0,250,50))  
