library(ggthemes)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(openxlsx)

#Load data
metadata_complete <- read.table("metadata_complete.txt", sep = "\t", quote = "#", header = T)
cazy <- read.table("CAZY_dbCAN2_node.txt", sep = "\t", quote = "#", header = T)
ent <- read.table("dirichlet_groups.txt", sep = "\t", quote = "#", row.names = 1)
matrix <- read.table("Count_matrix_all.txt", sep = "\t", quote = "#", header = T)
rownames(matrix) <- matrix$Gene
c_matrix <- matrix[matrix$Gene %in% cazy$Contig_ID,]
colnames(c_matrix)[1] <- "Contig_ID"
c_matrix <- c_matrix[order(c_matrix$Contig_ID),]

#Format CAZY data to analyze
cazy_m <- merge(cazy, c_matrix, by = "Contig_ID", all = T)
cazy_m[2113,c(1:10)]
write.table(cazy_m2, file = "Cazy_count_matrix.txt", sep = "\t")
cazy[cazy$Contig_ID %in%"NODE_1121_length_10222_cov_3.789401_1179_9",]
c_matrix[c_matrix$Contig_ID %in%"NODE_1121_length_10222_cov_3.789401_1179_9", c(1:6)]

cazy_m2 <- cazy_m[!cazy_m$CAZY %in% "", ]
cazy_m2 <- na.omit(cazy_m2)
cazy_m2 <- rowsum(cazy_m2[,c(3:ncol(cazy_m2))], as.character(cazy_m2$CAZY))
colnames(cazy_m2) <- sapply(colnames(cazy_m2), function(x){
  strsplit(x, "_")[[1]][1]
})

cazy_m2 <- cazy_m2[c(2:nrow(cazy_m2)),]
base_meta <- metadata_complete[metadata_complete$Time %in% "Baseline",]

base_meta[,1] <- sapply(1:nrow(base_meta), function(x){
  paste0("X", base_meta[x,1])
})
#CAZYs count matrix.
base_cazy <- cazy_m2[,base_meta$Life_ref]

write.table(base_cazy, "CAZY_count_matrix.txt", sep = "\t")
total_c <- apply(base_cazy, 2, sum)

base_cazy_rel <- as.data.frame(t(apply(base_cazy, 1, function(x){
  x*100/total_c
})))
mean_rel <- sort(apply(base_cazy_rel, 1, mean), decreasing = T)

##### Top 30 CAZY #### 
base_cazy_rel_30 <- base_cazy_rel[names(mean_rel[1:30]),]
df_base_cazy_rel_30 <- stack(base_cazy_rel_30)
colnames(df_base_cazy_rel_30) <- c("Values", "Samples")
df_base_cazy_rel_30$CAZy <- rownames(base_cazy_rel_30)
df_base_cazy_rel_30$CAZy_Group <- sapply(df_base_cazy_rel_30$CAZy, function(x){
  substr(x, 1, 2)
})
df_base_cazy_rel_30[df_base_cazy_rel_30$CAZy_Group %in% "CB","CAZy_Group" ] <- "CBM"
colnames(df_base_cazy_rel_30)
df_base_cazy_rel_30$CAZy <- factor(df_base_cazy_rel_30$CAZy, rev(names(mean_rel[1:30])))


#### CAZY GH29 GH 95 ###

GH_int <- base_cazy_rel[c("GH95", "GH29"),]

GH_int_t <- as.data.frame(t(GH_int))
rownames(GH_int_t)
base_meta$Life_ref
GH_int_t$Secretor <- base_meta[base_meta$Life_ref %in% rownames(GH_int_t),"HMtype2"]
GH_int_t$Samples <- rownames(GH_int_t)
GH_int_t <- na.omit(GH_int_t)
df_GH_int <- GH_int_t %>%
  melt(c("Secretor", "Samples"))
colnames(df_GH_int) <- c("Secretor", "Samples", "CAZY", "Values")
krust <- df_GH_int %>% 
  group_by(CAZY) %>%
  wilcox_test(Values ~ Secretor)
write.xlsx( krust, file = "KW-dunn_dbCAN2.xlsx", append = F, sheetName = "Krustal-Wallis")



##### CAZy GH #####
write.table(as.data.frame(t(mean_rel_GH)), "GH_mean_rel_abundance.txt", sep = "\t", row.names = F)
GH <- base_cazy_rel[grep("GH", rownames(base_cazy_rel)),]
mean_rel_GH <- mean_rel[names(mean_rel) %in% rownames(GH)]
GH_20 <- GH[names(mean_rel_GH)[1:20], ]
rownames(GH_20)
GH_20_t <- as.data.frame(t(GH_20))
GH_20_t$Community_type <- base_meta[base_meta$Life_ref %in% rownames(GH_20_t),"Community_type"]
GH_20_t$Samples <- rownames(GH_20_t)

df_GH_20 <- GH_20_t %>%
  melt(c("Community_type", "Samples"))
colnames(df_GH_20) <- c("Community_type", "Samples", "CAZY", "Values")

krust2 <- df_GH_20 %>% 
  group_by(CAZY) %>%
  kruskal_test(Values ~ Community_type) %>% 
  adjust_pvalue(method = "BH") %>%
  add_significance()
  
write.xlsx( krust2, file = "KW-dunn_dbCAN2.xlsx", append = F, sheetName = "Krustal-Wallis")

  stat.test <- df_GH_20 %>% 
    group_by(CAZY) %>%
    dunn_test(Values ~ Community_type, p.adjust.method = "BH")%>% 
    adjust_pvalue(method = "BH") %>%
    add_significance()

write.xlsx(stat.test, file = "KW-dunn_dbCAN2.xlsx", append = T, sheetName = "Dunn_test")
krust

df_GH_20_sig <- df_GH_20[as.character(df_GH_20$CAZY) %in% pull(krust[krust$p < 0.05,"CAZY"]),]
colnames(df_GH_20_sig)
df_GH_20_sig$Community_type <- factor(df_GH_20_sig$Community_type)
