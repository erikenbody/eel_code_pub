library(tidyverse)
library(data.table)
df.cov <- fread("data/indexcov/output/output-indexcov.ped", header = T)
df.cov2 <- fread("data/indexcov/output/output-indexcov.bed")
head(df.cov2)
df.covmeans <- colMeans(df.cov2[,-c(1:3)])

df.covsd <- apply(df.cov2[,-c(1:3)], 2, sd)

hist(df.covmeans,  breaks = 100)

hist(as.numeric(df.cov$mapped), breaks = 300)

df.cov %>% arrange(mapped) %>% head(n = 10) %>% select(sample_id)


#asrI15, atu10

# check outliers ----------------------------------------------------------

pca.outliers<-read.csv("output/pcangsd/pca_outliers.csv")
df.cov <- df.cov %>% 
  mutate(sample = gsub("Sample_TD-2544-", "",sample_id),
         outliers = ifelse(sample %in% pca.outliers$sample, "Y",NA))

df.cov$scaled.cov <- df.covmeans
df.cov$cov.sd <- df.covsd

df.cov %>% ggplot() + 
  geom_point(aes(x = unmapped, y = mapped, color = outliers), alpha = 0.5) +
  xlim(0,100000) + ylim(0,5e7)

df.cov %>% ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = outliers), alpha = 0.5) 

df.cov %>% ggplot() + 
  geom_point(aes(x = scaled.cov, y = cov.sd, color = outliers), alpha = 0.5) #+
  #xlim(0,1e8)

