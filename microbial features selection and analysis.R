## filter count #0
data_filter <- function(inputfile, num_S){
  data <- read.table(inputfile,sep = "\t", fill = T, header = T, row.names = 1)
  num_0 <- apply(data, 1, function(x) length(which(x == 0)))
  med <- apply(data, 1, mean)
  data_num <- cbind(data, as.data.frame(med), as.data.frame(num_0))
  colnames(data_num)[(ncol(data)+1):(ncol(data)+2)] <- c("median","count0")
  data_filter_0 <- data_num[which(data_num[,ncol(data_num)] <= as.numeric(num_S*0.8)),]
  data_filter <- data_filter_0[which(data_filter_0[,(ncol(data_filter_0)-1)] > 0.01),]
  return(data_filter[,1:(ncol(data))])
}
input_all <- "NC_species_C_filter_sparcc_per.txt"
data_filter_all <- data_filter(input_all, 124)

##PCA analysis
library(permute)
library(vegan)
library(devtools)
library(ggbiplot)
data_QS <- data_filter_all
group <- c(rep("Health",38),rep("NAFLD",86))
group_age <- c(as.vector(Health_meta[,'Age_group']),as.vector(NAFLD_meta[,'Age_group']))
group_age1 <- c(as.vector(Health_meta[,'AGE']),as.vector(NAFLD_meta[,'AGE']))
group_country <- c(rep("4",38),as.vector(NAFLD_meta[,'Race']))
group_sex <- c(as.vector(Health_meta[,'Gender']),as.vector(NAFLD_meta[,'Gender']))
data_log <- log2(100*data_QS+1)
decorana(t(data_log)) 
data_L <- t(data_log)
#data_L1 <- t(data_log1)
pca <- prcomp(data_L, scale = TRUE)
percent <- 100*pca$sdev^2/sum(pca$sdev^2)
label <- colnames(data_QS)
ggbiplot(pca,scale=TRUE, obs.scale = 1, var.scale = 1, ellipse = TRUE, groups = group_country,
         circle = TRUE, var.axes = FALSE, labels.size = 2, varname.size = 2,varname.abbrev = TRUE)
#distance.bray<-vegdist((data_L),method = 'bray')
#anosim.result<-anosim(distance.bray,group,permutations = 999)
##ANOSIM
anosim.result<-anosim(data_L,group_country,permutations = 9999,distance = "euclidean", strata = NULL,parallel = getOption("mc.cores"))
summary(anosim.result)
###alpha-diversity
library(vegan)
data_species <- data_filter_all
species_data <- t(data_species)
simp <- diversity(species_data, "simpson")
shan <- diversity(species_data, "shannon")
alpha <- data.frame(shan,simp)
write.table(alpha,"alpha_diversity_species_all.txt",sep = '\t',quote=TRUE,col.names= T, row.names = T)
## differental species
All_species <- data_filter_all
write.table(All_species,'NC_species_all.txt', sep = '\t',col.names= T, row.names = T)
#write.table(colnames(All_species)[39:124],"Case_sample.txt",sep = '\t')
data_NAFLD <- All_species[,which(substr(colnames(All_species),1,3)=="SRR")]
#data_Healt <- All_kegg[,c(which(substr(colnames(All_kegg),1,1)=="O"),which(substr(colnames(All_kegg),1,1)=="V"),which(substr(colnames(All_kegg),1,1)=="M"))]
data_Healt <- All_species[,c(which(substr(colnames(All_species),1,2)=="DE"),which(substr(colnames(All_species),1,2)=="FR"))]
table_mean <- cbind(apply(data_Healt,1,mean),apply(data_NAFLD,1,mean))
table_sd <- cbind(apply(data_Healt,1,sd),apply(data_NAFLD,1,sd))
table_med <- cbind(apply(data_Healt,1,median),apply(data_NAFLD,1,median))
logR <- apply(table_mean,1,function(x) log2(x[2]/x[1]))
table_test <- apply(All_species, 1, function(x) wilcox.test(x[1:ncol(data_Healt)], x[(ncol(data_Healt)+1):ncol(All_species)], paired = F)$p.value)
table_test <- cbind(table_test, p.adjust(table_test, method = "BH", n = length(table_test)))
table_diff <- table_test[table_test[,2]< 0.05,]
table_out <- cbind(All_species, table_mean, table_sd, logR, table_test)
colnames(table_out)[(ncol(All_species)+1):(ncol(All_species)+7)] <- c('H_mean','N_mean','H_sd','N_sd',
                                                                'lof2h/n','pvalue','FDR')
data_diff <- table_out[(table_out[,(ncol(table_out))] < 0.1),]
data_diff_log <- data_diff[abs(data_diff[,(ncol(data_diff)-2)]) > 1,]
data_interest <- data.frame(data_diff[,125:131])
data_interestl <- data.frame(data_diff_log[,125:131])
data_all <- data.frame(table_out[,125:131])
write.table(data_diff[,1:124],"NC_merge_species_diff.txt", sep = '\t',quote=TRUE,col.names= T, row.names = T)
write.table(data_interest,"NC_merge_species_diff_SD.txt", sep = '\t',quote=TRUE,col.names= T, row.names = T)
write.table(table_out[,125:131],"NC_merge_species_filter_SD.txt", sep = '\t',quote=TRUE,col.names= T, row.names = T)
##heatmap
library(pheatmap)
library("RColorBrewer")
data_heat <- read.table('NC_merge_species_diff_heatmap.txt',sep = "\t", fill = T, header = T, row.names = 1 )
data_heatlg <- log(100*data_heat[,1:2]+1)
colnames(data_heatlg) <- c("Health","NAFLD")
data_h <- cbind(data_heatlg,data_heat[,3:4])
write.table(data_h,"VIP_diff_heatmap.txt",sep = '\t',quote=TRUE,col.names= T, row.names = T)
annotation_row <- data.frame(data_heat[,3:4])
rownames(annotation_row) <- rownames(data_heatlg)
colnames(annotation_row) <- c('Class','Phylum')
anno_colors <- list(
  Phylum = c(Actinobacteria='#FFCC99',Bacteroidetes='#99CCFF',Firmicutes='#FF6666',Proteobacteria='#66CCCC'),
  Class = c(Actinobacteria='#FFCC99',Bacteroidia='#99CCCC',Bacilli='#FF99CC',Clostridia='#FFCCCC',Negativicutes='#FFFFCC',Erysipelotrichia='#CCCCFF',Betaproteobacteria='#CCFFFF',Deltaproteobacteria='#CCFFCC',
            Gammaproteobacteria='#009999')
)

pheatmap(data_heatlg,annotation_row = annotation_row,annotation_colors = anno_colors)
pheatmap(data_heatlg,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cluster_cols = FALSE,fontsize = 8,fontsize_row=11,cellwidth=35,cellheight = 10,
         annotation_row = annotation_row,annotation_colors = anno_colors,border_color = "white",scale = 'none')
anno_colors1 <- list(
  Phylum = c(Actinobacteria='#3E6E8A',Bacteroidetes='#76827E',Firmicutes='#FCE4B1',Proteobacteria='#E6A565'),
  Class = c(Actinobacteria='#001A4C',Bacteroidia='#2983CB',Bacilli='#E4E07F',Clostridia='#E8BD32',Negativicutes='#FFE782',Erysipelotrichia='#2F6083',Betaproteobacteria='#649AC3',Deltaproteobacteria='#58767C',
            Gammaproteobacteria='#A5CAD6')
)

mycolor <- colorRampPalette(c("#186D8C","#8DC8B2","#EADEB9","#F2AF66"))

pheatmap(data_heatlg,color = mycolor(100),cluster_cols = FALSE,fontsize = 12,fontsize_row=11,cellwidth=45,cellheight = 15.5,
         annotation_row = annotation_row,annotation_colors = anno_colors1,border_color = "grey60")
##bubble
data_VIP_diff_sta <- read.table("HN_species_VIP_diff_SD.txt",sep = "\t", fill = T, header = T, row.names = 1)
data_bu <- data.frame(data_VIP_diff_sta[,7:8])
colnames(data_bu) <- c('FDR','VIP')
logFDR=-log10(data_bu$FDR)
data<-cbind(data_bu,logFDR)
taxa <- data.frame(rownames(data))
data <- cbind(data, taxa)
x <- data.frame(rep(1,37))
colnames(x) <- 'level'
colnames(data)[4] <- 'taxa'
data <- cbind(data,x)
write.table(data,"bioplot.txt", sep = '\t',quote=TRUE,col.names= T, row.names = T)
data <- read.table("NC_merge_species_diff_Importance.txt", sep = '\t', row.names = 1, header=T)
idx <- order(data$level,decreasing = TRUE)
data$taxa <- factor(data$taxa, levels= data$taxa[idx]) 
ggplot(data,aes_string(x="axis",y="taxa", size='Importance', color='logFDR')) +
  geom_vline(xintercept = 1,linetype="dashed") +
  geom_point() + 
  scale_color_gradient(low="blue", high="red",limits=c(1,15))+
  scale_size_continuous( limits=c(0,8)) +
  theme(panel.background=element_blank())
scale_color_gradient(low="#186D8C", high="#A21E0E",limits=c(1,10))
##Phylum-ggallu
library(ggplot2)
library(ggalluvial)
library(reshape2)
library(tidyr)
data_per_filter <- function(inputfile){
  data <- read.table(inputfile,sep = "\t", fill = T, header = T, row.names = 1)
  data_sum <- apply(data,2,sum)
  data_per <- matrix(0,nrow(data),ncol(data))
  colnames(data_per) <- colnames(data)
  rownames(data_per) <- colnames(data)
  for (i in 1:nrow(data)){
    for (j in 1:ncol(data)){
      data_per[i,j] <- 100*data[i,j]/data_sum[j]
    }
  }
  num_0 <- apply(data_per, 1, function(x) length(which(x == 0)))
  med <- apply(data_per, 1, mean)
  data_num <- cbind(data_per, as.data.frame(med), as.data.frame(num_0))
  colnames(data_num)[(ncol(data)+1):(ncol(data)+2)] <- c("median","count0")
  data_filter_0 <- data_num[which(data_num[,ncol(data_num)] <= as.numeric(sample_num*0.8)),]
  data_filter <- data_filter_0[which(data_filter_0[,(ncol(data_filter_0)-1)] > 0.05),]
  return(data_filter[,1:(ncol(data))])
  return(data_per)
}
inputN <- "NAFLD_merge_abundance_phylum_C.txt"
data_per_N <- data_per_filter("NAFLD_merge_abundance_phylum_C.txt")
inputH <- "CRC_merge_abundance_phylum_C.txt"
data_per_N <- data_per_filter(inputH,38)
data_phylum <- read.table('NAFLD_merge_abundance_phylum_C.txt',sep = '\t', header = T, row.names = 1)
phylum_sum <- apply(data_phylum,2,sum)
data_phylum_per <- matrix(0,nrow(data_phylum),ncol(data_phylum))
colnames(data_phylum_per) <- colnames(data_phylum)
rownames(data_phylum_per) <- rownames(data_phylum)
for (i in 1:nrow(data_phylum)){
  for (j in 1:ncol(data_phylum)){
    data_phylum_per[i,j] <- 100*data_phylum[i,j]/phylum_sum[j]
  }
}
num_0 <- apply(data_phylum_per, 1, function(x) length(which(x == 0)))
med <- apply(data_phylum_per, 1, mean)
data_num <- cbind(data_phylum_per, as.data.frame(med), as.data.frame(num_0))
colnames(data_num)[(ncol(data_phylum)+1):(ncol(data_phylum)+2)] <- c("median","count0")
data_filter_0 <- data_num[which(data_num[,ncol(data_num)] <= as.numeric(86*0.8)),]
data_filter <- data_filter_0[which(data_filter_0[,(ncol(data_filter_0)-1)] > 0.05),]
data_phylum_H <- data_filter[,1:ncol(data_phylum)]
data_phylum_N <- data_filter[,1:ncol(data_phylum)]
data_phylum_H_mean <- data.frame(apply(data_phylum_H,1,mean))
data_phylum_N_mean <- data.frame(apply(data_phylum_N,1,mean))
colnames(data_phylum_H_mean) <-c('H_mean')
rownames(data_phylum_H_mean) <- rownames(data_phylum_H)
colnames(data_phylum_N_mean) <-c('N_mean')
rownames(data_phylum_N_mean) <- rownames(data_phylum_N)
write.table(data_phylum_N_mean,'Phylum_mean_N.txt',sep = '\t',quote=TRUE,col.names= T, row.names = T)
data_phy <- read.table("Phylum_mean.txt",sep = '\t', header = T, row.names = NULL)
melt_df = melt(data_phy)
ggplot(data = melt_df,
       aes(x = variable, weight = value, alluvium = Phylum)) +
  geom_alluvium(aes(fill = Phylum, colour = Phylum),
                alpha = .75,decreasing=TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text = element_text(size=12,colour = 'black')) +
  ggtitle("Phylum change among groups")
