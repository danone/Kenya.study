library(Maaslin2)
library(xlsx)

sp = read.table(file = "Metaphlan3/counts/Species_meta_counts.txt", header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
sp2[1:5, 1:5]
sp2<- rowsum(sp[,2:ncol(sp)], sp$Gene_ID)
colnames(sp2) <- sapply(colnames(sp2), function(x){
  strsplit(x, "_")[[1]][1]
})
df_input_data <- as.data.frame(t(sp2))
df_input_data[1:5, 1:5]
peak_microbiome = read.xlsx(file = "HMO_90 samples_28_02_2022.xlsx", sheetIndex = 2,  stringsAsFactors = FALSE)
rownames(peak_microbiome) <- peak_microbiome$Sample
peak_microbiome <- peak_microbiome[,-1]
peak_microbiome[is.na(peak_microbiome)] <- 0
colnames(peak_microbiome) <- gsub("HM", "", colnames(peak_microbiome))
complete_meta = read.table(file = "metadata_complete_211109.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df_input_metadata <- complete_meta[complete_meta$Time %in% "Baseline",]
rownames(df_input_metadata) <- paste0("X", rownames(df_input_metadata))

keep <- sapply(colnames(df_input_data), function(x){
  max(df_input_data[,x]) > 1000
})
df_input_data <- df_input_data[,keep]
log_peak_microbiome <-  as.data.frame(t(peak_microbiome))
colnames(log_peak_microbiome)<-gsub("'","",colnames(log_peak_microbiome))
colnames(log_peak_microbiome)<-gsub(" ","",colnames(log_peak_microbiome))
colnames(log_peak_microbiome)[1] <- "SL6"
colnames(log_peak_microbiome)[2] <- "SL3"
colnames(log_peak_microbiome)[5] <- "FL3"
colnames(log_peak_microbiome)[6] <- "FL2"
colnames(log_peak_microbiome)[7] <- "B4GL"
colnames(log_peak_microbiome)[23] <- "F3LNH"
colnames(log_peak_microbiome)[24] <- "F2LNH"
colnames(log_peak_microbiome)<-gsub("β","B",colnames(log_peak_microbiome))
colnames(log_peak_microbiome)<-gsub("-","",colnames(log_peak_microbiome))
colnames(log_peak_microbiome)
log_peak_microbiome<-log_peak_microbiome[order(rownames(log_peak_microbiome)),]

meta_hmo_base <- complete_meta[complete_meta$Time %in% "Baseline", ]
meta_hmo_base$Sample_ID <- sapply(as.character(meta_hmo_base$Sample_name), function(x){
  strsplit(x, " ")[[1]][2]
})
meta_hmo_base<-meta_hmo_base[order(meta_hmo_base$Sample_ID),]
meta_hmo_base <- meta_hmo_base[meta_hmo_base$Sample_ID %in% rownames(log_peak_microbiome), ]
sp_base <- df_input_data[rownames(df_input_data) %in% paste0("X", rownames(meta_hmo_base)),]
meta_hmo_base<-meta_hmo_base[order(rownames(meta_hmo_base)),]

meta_hmo_base[paste0("X", rownames(meta_hmo_base)) %in% rownames(sp_base), "Sample_ID"]
rownames(sp_base) <- meta_hmo_base[paste0("X", rownames(meta_hmo_base)) %in% rownames(sp_base), "Sample_ID"]

meta_hmo_base_peak <- meta_hmo_base
rownames(meta_hmo_base_peak) <- meta_hmo_base_peak$Sample_ID
meta_hmo_base_peak <- meta_hmo_base_peak[order(rownames(meta_hmo_base_peak)),]
rownames(meta_hmo_base_peak)
rownames(log_peak_microbiome)
meta_hmo_base_peak <- cbind(meta_hmo_base_peak, log_peak_microbiome)
colnames(meta_hmo_base_peak)
colnames(log_peak_microbiome)
k20 <- sort(apply(sp_base, 2, mean), decreasing = T)
k20 <- k20[c(1:20)]
sp_base_20 <- sp_base[,names(k20)]


sapply( colnames(log_peak_microbiome), function(x) {
  fit_data22 = Maaslin2(
    input_data = sp_base_20, 
    input_metadata = meta_hmo_base_peak, 
    output = paste0("peak_hmo_age_secretor", x),
    random_effects = c("age"),
    # normalization = "NONE",
    fixed_effects = c(x, "HMtype2"))
})

sapply( colnames(log_peak_microbiome), function(x) {
  fit_data22 = Maaslin2(
    input_data = sp_base_20, 
    input_metadata = meta_hmo_base_peak, 
    output = paste0("peak_hmo_age_", x),
    random_effects = c("age"),
    # normalization = "NONE",
    fixed_effects = c(x))
})
