library(RColorBrewer)
library(phyloseq)
library(wesanderson)

#### Panphlan loading #####


load("Panphlan_physeq.RData")
load("Panphlan_genes.RData")
load("Hmo_clusters_african.RData")

g_interest = read.table("HMO_genes_panphlan.txt", sep = "\t")
physeq2 <- physeq

SAM <- sample_data(physeq)

physeq = subset_samples(physeq, SAM != "Bifidobacterium_longum_subsp._infantis_CCUG_52486" & SAM != "Bifidobacterium_longum_subsp._infantis_157F" & SAM != "Bifidobacterium_longum_subsp._infantis_strain_NCTC13219_genome_assembly")

hmo_g = read.table("genes_hmo.txt", sep = "\t")
hmo_g <- hmo_g[order(hmo_g$V1),]

names <- dplyr::distinct(g_interest[,c(1,2)])
names <- names[order(names$V2),]
names$V2 <- as.character(names$V2)
names$V2[2] <- "araD1"


names2 <- data.frame( UniRef = names$V1, Name = names$V2, Gene = c(as.character(names[c(1:3),2] ), as.character(hmo_g$V3) ) )
names2 <- names2[order(names2$UniRef),]

OTU <- data.frame(otu_table(physeq)[rownames(otu_table(physeq)) %in% names2$UniRef,])
OTU<-OTU[order(rownames(OTU)),]


OTU$UniRef <- rownames(OTU)
t <- merge(names2, OTU, all.x = T )
t <- t[order(t$Gene), ]
# t<-t[-34,]
t<-t[-33,]

library(pheatmap)

#Prepare metadata

heat_col <- data.frame(SAM[,2, drop = F])
colnames(heat_col) <- "Genomes/Metagenomes"
heat_col[heat_col$`Genomes/Metagenomes` %in% "Bifidobacterium_longum_subsp._infatis_ATCC_15697","Genomes/Metagenomes"] <- "Bifidobacterium_longum_subsp._infantis"
rownames(heat_col) <- SAM$Sample

heat_col[heat_col$`Genomes/Metagenomes` %in% "African",] <- "Kenyan"
heat_col[heat_col$`Genomes/Metagenomes` %in% "Western",] <- "Swedish"

heat_row <- t[,3, drop =F]
heat_row$Genes <- c(rep("HMO_Cluster", 31), c("araA", "araD"))





color <-brewer.pal(n = 12, name = 'Paired')
p_color <- list(
  `Genomes/Metagenomes` = c(color[c(2,1,3,4 )]),
  Genes = c(color[c(7,8,9)])
)

names(p_color$`Genomes/Metagenomes`) <- unique(heat_col$`Genomes/Metagenomes`)
names(p_color$Genes) <- rev(unique(heat_row[,2]))



# Clustering
d <- dist(t(new_t[,c(4:ncol(new_t))]), method = "euclidean") # distance matrix
fit <- hclust(d, method="complete")
groups <- cutree(fit, k=4) 


g_all <- groups
g_all[g_all == "1"] <- "B. longum subsp. infantis"
g_all[g_all == "2"] <- "Putative B. longum subsp. infantis"
g_all[g_all == "3"] <- "B. longum subsp. longum"
g_all[g_all == "4"] <- "B. longum subsp. longum and B. longum subsp. infantis" 

heat_col2 <- heat_col
heat_col2[names(g_all),"Clusters"] <- g_all


heat_col2$`Genomes/Metagenomes` <- factor(heat_col2$`Genomes/Metagenomes`)
levels(heat_col2$`Genomes/Metagenomes`) <- c("B. longum subsp. infantis",
                                             "B. longum subsp. longum",
                                             "Kenyan", "Swedish")



color <-brewer.pal(n = 12, name = 'Paired')
color2 <- wes_palette("FantasticFox1", 5, type = c("discrete"))
p_color <- list(
  `Genomes/Metagenomes` = c(color[c(2,1,3,4 )]),
  Clusters = c(color2[c(2:5 )]),
  Genes = c(color[c(7,8,9)])
)

names(p_color$`Genomes/Metagenomes`) <- unique(heat_col2$`Genomes/Metagenomes`)
names(p_color$Clusters) <- unique(heat_col2$Clusters)[2:5]

names(p_color$Genes) <- rev(unique(heat_row[,2]))
p_color$Genes[1] <- '#FF7F00'


rnames <- new_t$Name
rnames[33] <- 'AraD'

heat_col2 <- heat_col2[row.names(heat_col2)%in%colnames(new_t),]

new_t2 <-new_t[order(new_t$Name),]

rnames <- new_t2$Name
rnames[2] <- 'AraD'

levels(heat_col2$`Genomes/Metagenomes`)[1] <- "B. longum subsp. infantis"
levels(heat_col2$`Genomes/Metagenomes`)[2] <- "B. longum subsp. longum"
heat_col2$Clusters <- factor(heat_col2$Clusters)
levels(heat_col2$Clusters)[1] <- "B. longum subsp. infantis"
levels(heat_col2$Clusters)[2] <- "B. longum subsp. longum" 
levels(heat_col2$Clusters)[3] <- "B. longum subsp. longum and B. longum subsp. infantis"
levels(heat_col2$Clusters)[4] <- "Putative B. longum subsp. infantis"
