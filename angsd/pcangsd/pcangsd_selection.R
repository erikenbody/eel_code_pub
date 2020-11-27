library(tidyverse)
library(patchwork)
library(RcppCNPy)

file.list <- list.files("data/pcangsd", full.names = TRUE, recursive = F)
file.list <- grep("selection.npy", file.list, value = TRUE)


chr.plots <- list()
df.list <- list()
for (chr in file.list){
  S <- npyLoad(chr) # Reads results from selection scan
  S.pos <- read.table(gsub(".selection.npy",".sites", chr))
  df.S <- as.data.frame(S)
  df.S <- cbind(df.S, S.pos)
  names(df.S) <- c("stat", "pos")
  df.S$row <- 1:nrow(df.S)
  df.S$chr <- gsub("eels_","",gsub(".selection.npy","",basename(chr)))
  df.S$pc1.pval <- 1 - pchisq(df.S$stat, 1)
  df.S$pc1.pval.rollmean <- zoo::rollmean(df.S$pc1.pval,50,fill=NA)
  #head(df.S, n = 50)
  chr.plots[[chr]]  <- ggplot(df.S) + 
    geom_point(aes(x = row, y = -log(pc1.pval))) +
    geom_line(aes(x = row, y = -log(pc1.pval.rollmean)), color = "red") +
    theme_bw() + ggtitle(unique(df.S$chr))
  
  df.list[[chr]] <- df.S
}


for (chr in file.list){
  print(chr.plots[[chr]])
  ggsave(paste0("output/pcangsd/", unique(df.S$chr), ".png"), width = 16, height = 6)
}




# define regions of interest ----------------------------------------------
names(df.list[["data/pcangsd/eels_SUPER_15.selection.npy"]]) <- c("stat", "pos", "row", "chr","pc1.pval", "pc1.pval.rollmean")
df.15 <- df.list[["data/pcangsd/eels_SUPER_15.selection.npy"]]

df.list[["data/pcangsd/eels_SUPER_15.selection.npy"]] %>% filter(row > 315000 & row < 340000) %>% 
  ggplot()+geom_point(aes(x = row, y = -log(pc1.pval))) + theme_bw() +
  geom_vline(aes(xintercept = 322600)) + geom_vline(aes(xintercept = 329300))
ggsave("output/pcangsd/super15_region_of_interest.png", width = 10, height = 8)


df.15 %>% 
  filter(row == 322600 | row == 329300) %>% select(pos)

#peak 2
df.list[["data/pcangsd/eels_SUPER_15.selection.npy"]] %>% filter(row > 755000 & row < 775000) %>% 
  ggplot()+geom_point(aes(x = row, y = -log(pc1.pval))) + theme_bw() +
  geom_vline(aes(xintercept = 757500)) + geom_vline(aes(xintercept = 770000))
ggsave("output/pcangsd/super15_region_of_interest_peak2.png", width = 10, height = 8)

ggsave("output/pcangsd/super15_region_of_interest_peak2.png", width = 10, height = 8)
df.15 <- df.list[["data/pcangsd/eels_SUPER_15.selection.npy"]]

df.15 %>% 
  filter(row == 757500 | row == 770000) %>% select(pos)


#peak 3
df.13 <- names(df.list[["data/pcangsd/eels_SUPER_13.selection.npy"]]) <- c("stat", "pos", "row", "chr","pc1.pval", "pc1.pval.rollmean")
df.13 <- df.list[["data/pcangsd/eels_SUPER_13.selection.npy"]]
df.13 %>% filter(row > 0 & row < 3.5e5) %>% 
  ggplot()+geom_point(aes(x = row, y = -log(pc1.pval))) + theme_bw() +
  geom_vline(aes(xintercept = 0)) + geom_vline(aes(xintercept = 3e5))
ggsave("output/pcangsd/super13_region_of_interest.png", width = 10, height = 8)

df.13 %>% 
  filter(row == 1 | row == 3e5) %>% select(pos)



# output ------------------------------------------------------------------

df.15 %>% write_csv("output/pcangsd/selection_SUPER_15.csv")
df.13 %>% write_csv("output/pcangsd/selection_SUPER_13.csv")


# manhattan ---------------------------------------------------------------
df.pval.genome <- bind_rows(df.list)
write_csv(df.pval.genome, "output/pcangsd/selection_stats_whole_genome_PC1.csv")
df.pval.genome %>% filter(-log(pc1.pval) > 1) %>% write_csv("output/pcangsd/selection_stats_whole_genome_PC1_logpca1pval_gr1.csv")
