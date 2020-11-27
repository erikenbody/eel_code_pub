library(tidyverse)
library(patchwork)

# -------------------------------------------------------------------------


file.list <- list.files("data/pcangsd", full.names = TRUE, recursive = F)
file.list <- grep(".cov", file.list, value = TRUE)

#outliers based on indexcov plots
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10")

df.list <- list()
plot.list <- list()
eval.list <- list()
for (pop in file.list){
  chr <- gsub(".cov","",basename(pop))
  pca.eel <- read.table(pop)
  pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted
  pca.names <- pca.names %>% filter(!grepl("ACA", pca.names$V1))

  pop.names <- read.csv("data/metadata/eel_pop_names.csv", stringsAsFactors = F)
  names(pca.names) <- "sample"
  pca.names$pop <- as.character(gsub('[[:digit:]]+', '', pca.names$sample))

  pca.names <- pca.names %>% left_join(pop.names, by = "pop")

  eigen.eel <- eigen(pca.eel, symmetric = T)
  eigenvectors.eel <- as.data.frame(eigen.eel$vectors)
  eigenvalues.eel <- eigen.eel$values

  eigenvectors.eel$sample <- pca.names$sample

  eigenvectors.eel$pop <- pca.names$pop
  eigenvectors.eel$Location <- pca.names$Location

  #drop outliers
  eigenvectors.eel <- eigenvectors.eel %>% filter(!sample %in% deviant.list)

  eigenvectors.eel[1:2,1:4]

  mycolors <- palette(c(
    rgb(64,124,109,  maxColorValue=255),
    rgb(215,158,190, maxColorValue=255),
    rgb(99,125,75, maxColorValue=255),
    rgb(118,91,136, maxColorValue=255),
    rgb(194,215,165, maxColorValue=255),
    rgb(62,116,153, maxColorValue=255),
    rgb(206,161,127, maxColorValue=255),
    rgb(149,179,224, maxColorValue=255),
    rgb(115,90,52, maxColorValue=255),
    rgb(118,203,198, maxColorValue=255),
    rgb(148,86,88, maxColorValue=255)))


  p1 <- ggplot() +
    geom_point(data = eigenvectors.eel, aes(x = V1, y = V2, color = Location)) +
    xlab(paste0("PC1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
    ylab(paste0("PC2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
    theme_bw() +
    scale_colour_manual(values = mycolors)+
    ggtitle(chr)

  plot.list[[pop]] <- p1
  eval.list[[pop]] <- eigenvalues.eel
  df.list[[pop]] <- eigenvectors.eel
}


wrap_plots(plot.list, ncol = 4)
ggsave("output/pcangsd/eel_all_pops_pca.png", height = 18, width = 18)


# plot 15 for figure ------------------------------------------------------
library(cowplot)
library(grid)
library(gridExtra)
df.super15 <- df.list[["data/pcangsd/eels_SUPER_15.cov" ]]
df.eval.super15 <- eval.list[["data/pcangsd/eels_SUPER_15.cov" ]]

p15 <- ggplot() +
  geom_point(data = df.super15, aes(x = V1, y = V2, color = Location), size = 2) +
  xlab(paste0("PC1", ": ", round(df.eval.super15[1],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(df.eval.super15[2], 1),"% variance")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text=element_text(size=21)) +
  scale_colour_manual(values = mycolors) +
  theme(legend.position = "none")
p15
ggsave("output/pcangsd/super_15_pca_plot.png", width = 8, height = 5, dpi = 320)

p.dummy <- p15 + theme(legend.position = c(0.5, 0.5))
legend <- cowplot::get_legend(p.dummy)


png("output/pcangsd/legened_locations.pdf", width = 8, height = 5, units = "in", res = 320)
grid.newpage()
grid.draw(legend)
dev.off()
