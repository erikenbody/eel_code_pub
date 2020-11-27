library(tidyverse)
library(patchwork)
library(cowplot)
library(grid)
library(gridExtra)
# -------------------------------------------------------------------------

#outliers based on indexcov plots
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10")


pca.eel <- read.table("data/pcangsd/subset_whole_genome/eels_all_subset_no_cananda.cov")

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

#eigenvectors.eel %>% select(V1, V2, sample) %>% View()
p1 <- ggplot() + 
  geom_point(data = eigenvectors.eel, aes(x = V1, y = V2, color = Location)) +
  xlab(paste0("PC1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
  ylab(paste0("PC2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
  theme_bw() +
  theme(text=element_text(size=21)) +
  theme(legend.position = "none") +
  scale_colour_manual(values = mycolors)# +
  #xlim(-0.1, 0.15) + ylim(-0.1, 0.2) #visualize without outliers
p1
ggsave("output/pcangsd/subset_whole_genome_pca_plot.pdf", width = 8, height = 5, dpi = 320, useDingbats = F)

# outliers ----------------------------------------------------------------
eigenvectors.eel %>% filter(V1 > 0.14) %>% select(sample) %>% write.csv("output/pcangsd/pca_outliers.csv")
# with canada -------------------------------------------------------------

pca.eel <- read.table("data/pcangsd/subset_whole_genome/eels_all_subset.cov")
pca.eel <- read.table("data/pcangsd/subset_whole_genome/eels_all_subset.cov")

pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 

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

#eigenvectors.eel %>% select(V1, V2, sample) %>% View()
p2 <- ggplot() + 
  geom_point(data = subset(eigenvectors.eel, Location!="Canada"), aes(x = V1, y = V2, color = Location)) +
  xlab(paste0("PC 1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
  ylab(paste0("PC 2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
  theme_bw() +
  theme(text=element_text(size=21)) +
  theme(legend.position = "none") +
  scale_colour_manual(values = mycolors)# +
#xlim(-0.1, 0.15) + ylim(-0.1, 0.2) #visualize without outliers
p2

ggsave("output/pcangsd/subset_whole_genome_pca_plot_canada_included_ZOOM.pdf", width = 8, height = 5, dpi = 320, useDingbats = F)

p3 <- ggplot() + 
  geom_point(data = eigenvectors.eel, aes(x = V1, y = V2, color = Location)) +
  xlab(paste0("PC 1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
  ylab(paste0("PC 2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
  theme_bw() +
  theme(text=element_text(size=21)) +
  theme(legend.position = "none") +
  scale_colour_manual(values = mycolors)# +
#xlim(-0.1, 0.15) + ylim(-0.1, 0.2) #visualize without outliers
p3

ggsave("output/pcangsd/subset_whole_genome_pca_plot_canada_included_FULL.pdf", width = 8, height = 5, dpi = 320, useDingbats = F)


p.dummy <- p2 + theme(legend.position = c(0.5, 0.5))
legend <- cowplot::get_legend(p.dummy)


png("output/pcangsd/legened_locations.pdf", width = 8, height = 5, units = "in", res = 320)
grid.newpage()
grid.draw(legend)
dev.off()



# flipped -----------------------------------------------------------------
#just add negative to x axis to flip
#eigenvectors.eel %>% select(V1, V2, sample) %>% View()
p5 <- ggplot() + 
  geom_point(data = subset(eigenvectors.eel, Location!="Canada"), aes(x = -V1, y = V2, color = Location)) +
  xlab(paste0("PC 1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
  ylab(paste0("PC 2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
  theme_bw() +
  theme(text=element_text(size=21)) +
  theme(legend.position = "none") +
  scale_colour_manual(values = mycolors)# +
#xlim(-0.1, 0.15) + ylim(-0.1, 0.2) #visualize without outliers
p5

ggsave("output/pcangsd/subset_whole_genome_pca_plot_canada_included_ZOOM_flip.pdf", width = 8, height = 5, dpi = 320, useDingbats = F)

p6 <- ggplot() + 
  geom_point(data = eigenvectors.eel, aes(x = -V1, y = V2, color = Location)) +
  xlab(paste0("PC 1", ": ", round(eigenvalues.eel[1],1),"% variance")) +
  ylab(paste0("PC 2", ": ", round(eigenvalues.eel[2], 1),"% variance")) +
  theme_bw() +
  theme(text=element_text(size=21)) +
  theme(legend.position = "none") +
  scale_colour_manual(values = mycolors)# +
#xlim(-0.1, 0.15) + ylim(-0.1, 0.2) #visualize without outliers
p6

ggsave("output/pcangsd/subset_whole_genome_pca_plot_canada_included_FULL_flip.pdf", width = 8, height = 5, dpi = 320, useDingbats = F)

