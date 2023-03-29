data.f2 <- data.frame(sample_data(physeq)[,c(11,14,15,19)])
colnames(data.f2) <- c('Age (Months)', 'pH', 'Calprotectin (ng/g)', 'Group')
# colnames(data.f2) <- c('Age', 'pH', 'Calprotectin', 'Group')

data.f2.melt <- reshape::melt(data.f2, id='Group')
data.f2.melt$value <- as.numeric(as.character(data.f2.melt$value))
colnames(data.f2.melt) <- c('Group', 'stat', 'value')
data.f2.melt$stat2 <- data.f2.melt$stat
levels(data.f2.melt$stat2)[3] <- expression(paste('Calprotectin (', mu, 'g/g)'))

#### Krutal walis ####

write.xlsx( data.f2.melt %>% 
              group_by(stat2) %>%
              kruskal_test(value ~ Group), file = "KW-dunn_pH-age.xlsx", append = F, sheetName = "Krustal-Wallis")

stat.test <- data.f2.melt %>% 
  group_by(stat2) %>%
  dunn_test(value ~ Group, p.adjust.method = "none")

write.xlsx(stat.test, file = "KW-dunn_pH-age.xlsx", append = T, sheetName = "Dunn_test")

stat.test
stat.test <- stat.test %>% add_y_position() 
stat.test$y.position[c(4,5)] <- c(8.5,8.9)
stat.test$p.adj.signif[c(4,5)] <- c("***",".")
stat.test


