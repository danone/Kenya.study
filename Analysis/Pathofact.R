library(phyloseq)
library(ggbiplot)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)

#Import data
meta = read.table(file='metadata_complete.txt', sep='\t', header = T, stringsAsFactors = F, fill = T)
meta[is.na(meta)] <- ""
meta <- meta[order(meta$Life_ref),]
base_meta <- meta[meta$Time %in% "Baseline",]

### African ARGs ###
#Upload count matrix
matrix <- read.table("../Count_matrix_all.txt", sep = "\t", quote = "#", header = T, stringsAsFactors = F)
c_m <- matrix
matrix <- c_m
colnames(matrix) <- c("Gene", sapply(2:ncol(matrix) , function(x){
  strsplit(colnames(matrix)[x], "_")[[1]][1]
}))
matrix<-matrix[order(matrix$Gene),]
rownames(matrix)<- matrix$Gene

###### matrix#####
## Create matrix
lengths <- read.table("../Counts_matrix_length.txt", sep = "\t", quote = "#", stringsAsFactors = F)
lengths <- lengths[order(lengths$V1),]

matrix <- matrix[,order(colnames(matrix))]


## Uploading the ARG information
patho <- read.table("../PathoFact_ORFs_clusters_nt_predictions.tsv", sep = "\t", quote = "#", header = T, stringsAsFactors = F)
patho<-patho %>% distinct(Contig, .keep_all = TRUE)
patho<-patho[order(patho$Contig), ]

####### Working with the ARG values

ARG <- patho[patho$ARG != "-",]
length(unique(patho$ARG))
ARG_c <- matrix[matrix$Gene %in% ARG$Contig,]
ARG_c <- ARG_c[order(ARG_c$Gene),]
ARG_sample<- apply(ARG_c[,c(2:ncol(ARG_c))], 2, function(x){
  length(x[x > 1])
  
})

ARG_sample<-ARG_sample[order(names(ARG_sample))]
ARG_sample<-ARG_sample[names(ARG_sample)%in% paste0("X", meta[meta$Time %in% "Baseline", "Life_ref"])]

base_matrix <- matrix[,colnames(matrix) %in% paste0("X", meta[meta$Time %in% "Baseline", "Life_ref"])]
total_counts <- apply(base_matrix, 2, sum)
ARG_sample<-ARG_sample[order(names(ARG_sample))]
ARG_sample_rel <- ARG_sample *100/total_counts
sp_abun_rel_df <- stack(as.data.frame(t(ARG_sample_rel)))
sp_abun_rel_df$Enterotype <- meta[paste0("X", meta$Life_ref) %in% names(ARG_sample_rel), "entero"]
sp_abun_rel_df$Pathofact <- "ARG"
sp_abun_rel_df$ind
meta[paste0("X", meta$Life_ref) %in% names(ARG_sample_rel), "Life_ref"]




#### TOxins #####



contig <- patho[,c("ORF_ID", "Toxin_confidence_level", "Contig")]
rownames(contig) <- patho$ORF

toxin <- read.table("../Toxin_gene_library_ORFs_clusters_nt_report.tsv", sep = "\t", quote = "#", header = T, stringsAsFactors = F)
toxin<-toxin %>% distinct(ORF, .keep_all = TRUE)
rownames(toxin) <- toxin$ORF
contig<-contig[contig$ORF_ID %in% toxin$ORF_ID,]
contig<-contig[order(contig$ORF_ID),]
toxin <- toxin[toxin$ORF_ID %in% contig$ORF_ID,]
toxin<-toxin[order(toxin$ORF_ID),]
toxin$Contig <- contig$Contig
NAME_c <- matrix[matrix$Gene %in% toxin$Contig,]
NAME_c <- NAME_c[order(NAME_c$Gene),]
NAME_c_all <- NAME_c
NAME_c_all$ID <- contig[contig$Contig %in% NAME_c$Gene, "Toxin_confidence_level"]
NAME_c_all <- NAME_c_all

###### Obtain Total toxins per sample ########

toxin_sample<- apply(NAME_c_all[,c(2:(ncol(NAME_c_all)-1))], 2, function(x){
  length(x[x > 1])
  
})

length(toxin_sample)
toxin_sample<-toxin_sample[names(toxin_sample)%in% paste0("X", meta[meta$Time %in% "Baseline", "Life_ref"])]
toxin_sample<-toxin_sample[order(names(toxin_sample))]

toxin_sample_rel <- toxin_sample *100/total_counts
toxin_abun_rel_df <- stack(as.data.frame(t(toxin_sample_rel)))
toxin_abun_rel_df$Enterotype <- meta[paste0("X", meta$Life_ref) %in% names(ARG_sample), "entero"]
toxin_abun_rel_df$Pathofact <- "Toxin"


##### VF #### 


VF <- read.table("../Virulence_prediction_ORFs_clusters_nt_report.tsv", sep = "\t", quote = "#", header = T, stringsAsFactors = F)
VF <- VF[!VF$Virulence_confidence_level %in% "-",]

contig_VF <- patho[,c("ORF_ID", "Virulence_confidence_level", "Contig")]
contig_VF<-contig_VF[contig_VF$ORF_ID %in% VF$ORF_ID,]
contig_VF<-contig_VF[order(contig_VF$ORF_ID),]

VF <- VF[VF$ORF_ID %in% contig_VF$ORF_ID,]
VF<-VF[order(VF$ORF_ID),]
VF$Contig <- contig_VF$Contig
VF_c <- matrix[matrix$Gene %in% VF$Contig,]
VF_c <- VF_c[order(VF_c$Gene),]
VF_c_all <- VF_c
VF_c_all$ID <- contig_VF[contig_VF$Contig %in% VF_c$Gene, "Toxin_confidence_level"]

VF_sample<- apply(VF_c_all[,c(2:(ncol(VF_c_all)))], 2, function(x){
  length(x[x > 1])
  
})

VF_sample<-VF_sample[names(VF_sample)%in% paste0("X", meta[meta$Time %in% "Baseline", "Life_ref"])]
VF_sample<-VF_sample[order(names(VF_sample))]
VF_sample_rel <- VF_sample *100/total_counts
vf_abun_rel_df <- stack(as.data.frame(t(VF_sample_rel)))
vf_abun_rel_df$Enterotype <- meta[paste0("X", meta$Life_ref) %in% names(ARG_sample), "entero"]
vf_abun_rel_df$Pathofact <- "VF"

all <- rbind(vf_abun_rel_df, toxin_abun_rel_df)
all2 <- rbind(all, sp_abun_rel_df)


#### Krutal walis ####
library(rstatix)
write.xlsx( all2 %>% 
  group_by(Pathofact) %>%
  kruskal_test(values ~ Enterotype), file = "KW-dunn_PathoFact.xlsx", append = F, sheetName = "Krustal-Wallis")

stat.test <- all2 %>% 
  group_by(Pathofact) %>%
  dunn_test(values ~ Enterotype, p.adjust.method = "none")

