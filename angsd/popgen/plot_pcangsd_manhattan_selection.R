library(tidyverse)
library(patchwork)
library(RcppCNPy)

df.pval.genome <- read.csv("output/pcangsd/selection_stats_whole_genome_PC1_logpca1pval_gr1.csv")
head(df.pval.genome)
df.pval.genome$chr_labels <- as.numeric(gsub("SUPER_","",df.pval.genome$chr))
chr_order <- df.pval.genome %>% group_by(chr) %>% 
  summarise(chr_labels = max(chr_labels)) %>% arrange(chr_labels)

df.pval.genome$chr_ordered <- factor(df.pval.genome$chr, levels = chr_order$chr)

df.in <- df.pval.genome %>% arrange(chr_labels)
df.in$row <- 1:nrow(df.in)
chr_breaks <- df.in %>% filter(!is.na(row)) %>% 
  #mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
  group_by(chr_ordered, chr_labels) %>% 
  dplyr::summarise(chr_breaks = mean(row))

chrom.colors <- data.frame(chr_ordered = grep("SUPER", unique(df.in$chr_ordered), value = T),
                           color.num = rep(1:2,length(grep("SUPER", unique(df.in$chr_ordered))))) %>% 
  distinct(chr_ordered, .keep_all = T)

df.in2 <- df.in %>% #mutate(row = 1:n()) %>% 
  left_join(chrom.colors, by = "chr_ordered") %>% 
  mutate(color.num = as.factor(color.num))

df.in2 %>% filter(chr_labels != "unknown" & !is.na(row)) %>% #filter(chr_ordered == "SUPER_19") %>% 
  ggplot(aes(x = row, y = -log(pc1.pval), col = color.num)) + theme_bw() +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, color = "black"),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2)) +
  geom_point(size=0.9,shape=20,stroke=0.2) +
  scale_color_manual(values=rep(c("grey30","grey70"))) +
  #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
  scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                     labels = function(labels) {
                       sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                     }) +
  ylab(NULL) 
ggsave("output/pcangsd/manhattan_pcangsd.png", width = 12, height = 2, dpi = "retina")


# zoom in -----------------------------------------------------------------

df.in2 %>% filter(chr == "SUPER_15") %>% 
  separate(pos, into = c(NA,NA,"position")) %>% 
  ggplot(aes(x = position, y = -log(pc1.pval))) + theme_bw() +
  geom_line(aes(x = position, y = -log(pc1.pval.rollmean)), color = "red") +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, color = "black"),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2)) +
  geom_point(size=0.9,shape=20,stroke=0.2) +
  labs(y= "-log(PC1 pvalue)") 
