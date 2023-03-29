peaks <- read.table('Copie de DataAnalysisTable_rqGFP_145_Kenyan_HM_samples_20210618(32761)_2.csv', header = F, sep = ';', dec=',', quote = '', fill = T)

peak_intensity <- peaks[,c(1,2,seq(3,1307, by = 9))]
peak_intensity_noisenorm <- peaks[,c(1,2,seq(4,1307, by = 9))]
peak_area <- peaks[,c(1,2,seq(5,1307, by = 9))]
peak_intensity_ISnorm <- peaks[,c(1,2,seq(6,1307, by = 9))]
peak_area_ISnorm <- peaks[,c(1,2,seq(7,1307, by = 9))]
peak_intensity_ISnorm_LOQ <- peaks[,c(1,2,seq(8,1307, by = 9))]
peak_area_ISnorm_LOQ <- peaks[,c(1,2,seq(9,1307, by = 9))]
peak_intensity_ISnorm_LOQ_sqrt <- peaks[,c(1,2,seq(10,1307, by = 9))]
peak_area_ISnorm_LOQ_sqrt <- peaks[,c(1,2,seq(11,1307, by = 9))]#Use this one

keep = peaks$V2!=''&peaks$V2!='OS-Ladder-associated'
restrict <- c('HM002', 'HM008', 'HM036', 'HM235', 'HM249','HM260')

names = read.table('/vol1/P_192543_6/9_statistics/6_HMO_peaks/names.txt', header = F, sep='\t', stringsAsFactors = F)
names_LOQ = read.table('/vol1/P_192543_6/9_statistics/6_HMO_peaks/names_LOQ.txt', header = F, sep='\t', stringsAsFactors = F)

colnam <- c('PeakID','Peak',as.character(names[1,]))
annot <- as.factor(as.character(names[2,]))

colnames(peak_intensity) <- colnam
colnames(peak_intensity_noisenorm) <- colnam
colnames(peak_area) <- colnam
colnames(peak_intensity_ISnorm) <- colnam
colnames(peak_area_ISnorm) <- colnam

colnames(peak_area_ISnorm_LOQ_sqrt) <- colnam


peak_area_ISnorm_LOQ_sqrt_2 <- (peak_area_ISnorm_LOQ_sqrt[keep,-c(colnam%in%restrict)])
row.names(peak_area_ISnorm_LOQ_sqrt_2) <- peak_area_ISnorm_LOQ_sqrt_2$Peak
peak_area_ISnorm_LOQ_sqrt_2 <- peak_area_ISnorm_LOQ_sqrt_2[,-1]
restrictnames <- colnames(peak_area_ISnorm_LOQ_sqrt_2)%in%restrict
peak_area_ISnorm_LOQ_sqrt_2 <- peak_area_ISnorm_LOQ_sqrt_2[,c(!restrictnames)]

annot <- annot[c(!restrictnames)]

annot2 <- annot 
levels(annot2) <- c('Secretor', 'Non secretor', 'Secretor', 'Non secretor')


plot(as.numeric(peak_area_ISnorm_LOQ_sqrt_2[1,])~annot)

write.table(file = '/vol1/P_192543_6/9_statistics/6_HMO_peaks/mean_peaks_LOQ_area.txt', 
            data.frame(
  TypeI=apply(peak_area_ISnorm_LOQ_sqrt_2[,annot=='HM Type I'], 1, mean, na.rm =T),
  TypeII=apply(peak_area_ISnorm_LOQ_sqrt_2[,annot=='HM Type II'], 1, mean, na.rm =T),
  TypeIII=apply(peak_area_ISnorm_LOQ_sqrt_2[,annot=='HM type III'], 1, mean, na.rm =T),
  TypeIV=apply(peak_area_ISnorm_LOQ_sqrt_2[,annot=='HM Type IV'], 1, mean, na.rm =T)
  ), sep='\t', quote=F, col.names = T, row.names = T)



plot(as.numeric(peak_area_ISnorm_LOQ_sqrt_2['3-FL',])~annot)


library(reshape2)


peak_area_ISnorm_LOQ_sqrt_2$name <- row.names(peak_area_ISnorm_LOQ_sqrt_2)
aqm <- melt(peak_area_ISnorm_LOQ_sqrt_2, id=c('name'), na.rm=T)
aqm$value <- as.numeric(aqm$value)
annot.df <- data.frame(annot)
row.names(annot.df) <- colnames(peak_area_ISnorm_LOQ_sqrt_2)[1:139]
aqm$annot <- annot.df[as.character(aqm$variable),1]

library(ggplot2)

aqm_IS <- aqm[aqm$name!='Internal Standard (IS)', ]
aqm_original <- aqm
aqm <- aqm_IS
levels(aqm$annot) <- c('HM-I', 'HM-II', 'HM-III', 'HM-IV')
secret <- aqm$annot
aqm$name[aqm$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
aqm$name[aqm$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
nameHMO_sort <- names(sort(apply(peak_area_ISnorm_LOQ_sqrt_2[,1:139],1, function(x){sum(as.numeric(x), na.rm = T)}), decreasing = T))
nameHMO_sort <- nameHMO_sort[-4]
nameHMO_sort[c(4,7)] <- c('LNFP II', 'LNFP I')

aqm$name <- factor(aqm$name, levels = nameHMO_sort)
levels(secret) <- c('Secretor', 'Non secretor', 'Secretor', 'Non secretor')
pdf('peacks_facet_noIS_LOQ_area.pdf', width = 11)
ggplot(aqm, aes(x=annot, y=value, fill=secret)) + geom_boxplot()  +scale_fill_manual(values = c('cornflowerblue', 'orchid1')) +
  facet_wrap( ~ name, ncol = 5, scales = 'free_y')
dev.off()
svg('peacks_facet_noIS_LOQ_area.svg', width = 11)
ggplot(aqm, aes(x=annot, y=value, fill=secret)) + geom_boxplot()  + scale_fill_manual(values = c('cornflowerblue', 'orchid1')) +
  facet_wrap( ~ name, ncol = 5, scales = 'free_y')
dev.off()

levels(secret) <- c('Secretor', 'Non secretor', 'Secretor', 'Non secretor')
pdf('peacks_facet_2_noIS_LOQ_area.pdf', width = 11)
ggplot(aqm, aes(x=secret, y=value, fill=secret)) + geom_boxplot()  +scale_fill_manual(values = c('cornflowerblue', 'orchid1')) +
  facet_wrap( ~ name, ncol = 5, scales = 'free_y')  
dev.off()
svg('peacks_facet_2_noIS_LOQ_area.svg', width = 11)
ggplot(aqm, aes(x=secret, y=value, fill=secret)) + geom_boxplot()  + scale_fill_manual(values = c('cornflowerblue', 'orchid1')) +
  facet_wrap( ~ name, ncol = 5, scales = 'free_y')
dev.off()

library(xlsx)

write.xlsx(
t(apply(!apply(peak_area_ISnorm_LOQ_sqrt_2, 1, is.na), 2, function(x) {c(Secretor=sum(x[annot2=='Secretor']), NonSecretor=round(sum(x[annot2=='Non secretor'])/39*100))})),
sheetName = 'Secretor group', file = 'Secretion_detection_percent.xlsx')

write.xlsx(
t(apply(!apply(peak_area_ISnorm_LOQ_sqrt_2, 1, is.na), 2, function(x) {c(TypeI=round(sum(x[annot=='HM Type I'])/73*100),
                                                                     TypeII=round(sum(x[annot=='HM Type II'])/31*100),
                                                                     TypeIII=round(sum(x[annot=='HM type III'])/27*100),
                                                                     TypeIV=round(sum(x[annot=='HM Type IV'])/8*100)
                                                                     )})),
sheetName = 'Secretor type', file = 'Secretion_detection_percent.xlsx', append = T)

aqm_2 <- aqm[order(aqm$annot),]


names_1 <- (apply(peak_area_ISnorm_LOQ_sqrt_2[,1:139],2, sum, na.rm=T))
names_1_1 <- names_1[annot=='HM Type I']
names_1_1 <- sort(names_1_1)
names_1_2 <- names_1[annot=='HM Type II']
names_1_2 <- sort(names_1_2)
names_1_3 <- names_1[annot=='HM type III']
names_1_3 <- sort(names_1_3)
names_1_4 <- names_1[annot=='HM Type IV']
names_1_4 <- sort(names_1_4)
names_sorted <- c(names(names_1_1), names(names_1_2),names(names_1_3),names(names_1_4))


aqm_2$variable <- factor(aqm_2$variable, levels=names_sorted)

aqm_2$name[aqm_2$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
aqm_2$name[aqm_2$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'

aqm_2 <- aqm_2[(aqm_2$name!="Internal Standard (IS)"),] 


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

pdf('peacks_barplot_2_area_LQD.pdf', width = 10)
ggplot(aqm_2, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
  scale_fill_manual(values = mycolors) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

svg('peacks_barplot_2_area_LQD.svg', width = 12)
ggplot(aqm_2, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
  scale_fill_manual(values = mycolors) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


peak_area_ISnorm_LOQ_sqrt_4 <- peak_area_ISnorm_LOQ_sqrt_2[,1:139]
databarplot0 <- data.frame(Secretor=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM Type I'|annot=='HM type III'], na.rm=T))}),
           NonSecretor=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM Type II'|annot=='HM Type IV'], na.rm=T))}),
           name=row.names(peak_area_ISnorm_LOQ_sqrt_4)
)

databarplot1 <- data.frame(Secretor=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM Type I'|annot=='HM type III'], na.rm=T))}),
                           NonSecretor=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM Type II'|annot=='HM Type IV'], na.rm=T))}),
                           name=row.names(peak_area_ISnorm_LOQ_sqrt_4)
)
           
 aqm_3 <- melt(databarplot0, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
 aqm_3 <- aqm_3[(aqm_3$name!="Internal Standard (IS)"),] 
 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 levels(aqm_3$name)[2] <- "2'-FL"
 levels(aqm_3$name)[24] <- "β4-GL"
 
 HMO_barplot_2groups <-  ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors) +  theme_pubr(legend='right') +ylab('mean of IS-normalized Peak Area') +xlab('')+
   guides(fill=guide_legend(title=""))
 
 svg('barplot_secretor_2_mean_area_LQD_goodcolors_goodnames.svg')
 HMO_barplot_2groups
 dev.off()
 
 pdf('barplot_secretor_2_mean_area_LQD_goodcolors_goodnames.pdf')
 HMO_barplot_2groups
 dev.off()
 ##
 svg('barplot_secretor_2_mean_area_LQD.svg')
ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
dev.off()

save(HMO_barplot_2groups, file = 'HMO_barplot_2groups.RData')

 aqm_3 <- melt(databarplot1, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
 aqm_3 <- aqm_3[(aqm_3$name!="Internal Standard (IS)"),] 
 
 svg('barplot_secretor_2_sum_area_LQD.svg')

 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 svg('barplot_secretor_2_sum_area_LQD_goodcolors.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 ##
 
 databarplot3 <- data.frame(HM_I=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM Type I'], na.rm=T))}),
                            HM_II=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM Type II'], na.rm=T))}),
                            HM_III=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM type III'], na.rm=T))}),
                            HM_IV=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(mean(x[annot=='HM Type IV'], na.rm=T))}),
                            name=row.names(peak_area_ISnorm_LOQ_sqrt_4)
 )
 
 databarplot4 <- data.frame(HM_I=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM Type I'], na.rm=T))}),
                            HM_II=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM Type II'], na.rm=T))}),
                            HM_III=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM type III'], na.rm=T))}),
                            HM_IV=apply(peak_area_ISnorm_LOQ_sqrt_4, 1, function(x){(sum(x[annot=='HM Type IV'], na.rm=T))}),
                            name=row.names(peak_area_ISnorm_LOQ_sqrt_4)
 )
 library(reshape)
 aqm_3 <- melt(databarplot3, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
 aqm_3 <- aqm_3[(aqm_3$name!="Internal Standard (IS)"),] 
 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 svg('barplot_secretor_3_mean_area_LQD_goodcolors.svg')
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 ##
 
 
 svg('barplot_secretor_3_mean_area_LQD.svg')
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 aqm_3 <- melt(databarplot4, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
 aqm_3 <- aqm_3[(aqm_3$name!="Internal Standard (IS)"),] 
 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 svg('barplot_secretor_3_sum_area_LQD_goodcolors.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 ##
 
 
 svg('barplot_secretor_3_sum_area_LQD.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors)
 dev.off()
 
 
 databarplot5 <- databarplot3
 databarplot5 <- databarplot5[row.names(databarplot5)!=("Internal Standard (IS)"),]
 databarplot5[is.na(databarplot5)] <- 0
 databarplot5[,1:4] <-  t(t(databarplot5[,1:4])/(apply( databarplot5[,1:4], 2, sum))*100)
 
 aqm_3 <- melt(databarplot5, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'

 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 svg('barplot_secretor_3_relative_LOQ_area_LQD_goodcolors.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors) + ylab('Value(%)')
 dev.off()
 ##
 svg('barplot_secretor_3_relative_LOQ_area_LQD.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors) + ylab('Value(%)')
 dev.off()
 
 
 
 databarplot5 <- databarplot1
 databarplot5 <- databarplot5[row.names(databarplot5)!=("Internal Standard (IS)"),]
 databarplot5[is.na(databarplot5)] <- 0
 databarplot5[,1:2] <-  t(t(databarplot5[,1:2])/(apply( databarplot5[,1:2], 2, sum))*100)
 
 aqm_3 <- melt(databarplot5, id=c('name'), na.rm=F)
 aqm_3$name <- as.character( aqm_3$name)
 
 aqm_3$name[aqm_3$name=="LNFP I (+ minor and currently unknown peak)"] <- 'LNFP I'
 aqm_3$name[aqm_3$name=="LNFP II (+ minor and currently unknown peak)"] <- 'LNFP II'
 

 aqm_3$name <- factor(aqm_3$name, levels=levels(aqm_2$name))
 svg('barplot_secretor_4_relative_LOQ_area_LQD_goodcolors.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors) + ylab('Value(%)')
 dev.off()
 ##
 
 svg('barplot_secretor_4_relative_LOQ_area_LQD.svg')
 
 ggplot(aqm_3, aes(x=variable, y=value, fill=name)) + geom_bar(stat='identity')+
   scale_fill_manual(values = mycolors) + ylab('Value(%)')
 dev.off()
