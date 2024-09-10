data <- read.delim('803sample_info.txt',header = T)

con2cou <- unique(data[,c("continent","Country")])
con2cou_num <- aggregate(con2cou$Country,list(con2cou$continent),length)
con2cou_num

con2spe <- unique(data[,c("continent","Species_Name")])
con2spe_num <- aggregate(con2spe$Species_Name,list(con2spe$continent),length)
con2spe_num

con2sam <- unique(data[,c("continent","fq_name")])
con2sam_num <- aggregate(con2sam$fq_name,list(con2sam$continent),length)
con2sam_num

con2w_c <- unique(data[,c("continent","Wild.Captive","fq_name")])
con2w_c$habitat_new <- ifelse(con2w_c$Wild.Captive == "Wild","Wild","Captive")
con2w_c_num <- aggregate(con2w_c$fq_name,list(con2w_c$continent,con2w_c$habitat_new),length)
con2w_c_num

data <- read.delim('39622genome_info.txt',header = T)
con2gen <- data[,c("continent","genome")]
con2gen_num <- aggregate(con2gen$genome,list(con2gen$continent),length)
con2gen_num
