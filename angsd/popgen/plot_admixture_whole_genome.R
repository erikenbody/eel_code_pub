library(tidyverse)
library(patchwork)
library(RcppCNPy)
library(ggthemes)
library(forcats)

location.order <- c("Canada", "Ireland", "England", "France", "Portugal", "Morocco", "Tunisia",
                    "AST", "SMS", "ASRI", "Lithuania")
#pop.names
#methods are largely indistinguishable
# pcangsd method ----------------------------------------------------------

S.admix <- npyLoad("data/pcangsd/subset_whole_genome/eels_all_subset.admix.Q.npy") # Reads results from selection scan
pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10") #low coverage

df.admix <- as.data.frame(S.admix)
df.admix$sampleID <- pca.names$V1

df.admix <- df.admix %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -sampleID) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.admix$loc <- factor(ifelse(grepl("ACA", df.admix$sampleID), "American", "European"))
df.admix$pop <- as.character(gsub('[[:digit:]]+', '', df.admix$sampleID))
pop.names <- read.csv("data/metadata/eel_pop_names.csv", stringsAsFactors = F)
df.admix <- df.admix %>% left_join(pop.names, by = "pop")

df.admix <- df.admix %>% filter(!sampleID %in% deviant.list)
#https://luisdva.github.io/rstats/model-cluster-plots/

ggplot(df.admix, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf"), guide = F) +
  scale_color_manual(values = c("#ef8a62", "#67a9cf"), guide = F)

ggsave("output/pcangsd/admixture_with_canada_whole_genome.png", width = 10, height = 4)
ggsave("output/pcangsd/admixture_with_canada_whole_genome.pdf", width = 10, height = 4, dpi = 320)


ggplot(df.admix, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Location), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  #scale_fill_gdocs(guide = FALSE) +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf"), guide = F)
ggsave("output/pcangsd/admixture_with_canada_whole_genome_location.png", width = 10, height = 4)
ggsave("output/pcangsd/admixture_with_canada_whole_genome_location.pdf", width = 10, height = 4, dpi = 320)


# NGSadmix method ---------------------------------------------------------
#files say super 1 but are whole genome

S.admix <- read.table("data/ngsadmix/subset_whole_genome/ngs_admix_k2_super1.qopt") # Reads results from selection scan
pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10")

df.admix.k2 <- as.data.frame(S.admix)
df.admix.k2$sampleID <- pca.names$V1

df.admix.k2 <- df.admix.k2 %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -sampleID) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.admix.k2$loc <- factor(ifelse(grepl("ACA", df.admix.k2$sampleID), "American", "European"))
df.admix.k2$pop <- as.character(gsub('[[:digit:]]+', '', df.admix.k2$sampleID))
pop.names <- read.csv("data/metadata/eel_pop_names.csv", stringsAsFactors = F)
df.admix.k2 <- df.admix.k2 %>% left_join(pop.names, by = "pop")

df.admix.k2 <- df.admix.k2 %>% filter(!sampleID %in% deviant.list)
#https://luisdva.github.io/rstats/model-cluster-plots/
df.admix.k2$Location <- factor(df.admix.k2$Location, levels = location.order)


p.ngsadmix.k2 <- ggplot(df.admix.k2, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  #geom_col(color = "gray", size = 0.1) +
  geom_col(aes(color = popGroup), size = 0.1) +
  facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("grey37", "#ef8a62"), guide = F) +
  scale_color_manual(values = c("grey27", "#ef8a62"), guide = F)
p.ngsadmix.k2

ggsave("output/ngsadmix/admixture_with_canada_whole_genome_k2.png", width = 10, height = 4)
ggsave("output/ngsadmix/admixture_with_canada_whole_genome_k2.pdf", width = 10, height = 4)

p.ngsadmix.k2.loc <- ggplot(df.admix.k2, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  #theme(strip.text.x = element_text(angle = 45)) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("#67a9cf", "#ef8a62"), guide = F) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62"), guide = F)


# plot proportion euro ----------------------------------------------------

american.ancestry <- df.admix.k2 %>% filter(loc == "European" & popGroup == 1) %>%#pop group 1 in this case is American
  group_by(Location) %>% 
  summarise(american_ancestry = 100 * round(mean(prob), 3))

p.minor.ancestry <- df.admix.k2 %>% filter(loc == "European" & popGroup == 1 | loc == "American" & popGroup == 2) %>% #View()
  ggplot() + geom_boxplot(aes(x = Location, y = prob, color = loc)) +
  scale_color_manual(values = c("#ef8a62", "#67a9cf"), guide = F) +
  theme_minimal() + 
  labs(x = NULL, title = "K=2", y = "Minor ancestry proportion") 


df.admix.k2 %>% filter(loc == "European" & popGroup == 1) %>% arrange(-prob)

# k=3 ---------------------------------------------------------
#files say super 1 but are whole genome

S.admix <- read.table("data/ngsadmix/subset_whole_genome/ngs_admix_k3_super1.qopt") # Reads results from selection scan
pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10")

df.admix.k3 <- as.data.frame(S.admix)
df.admix.k3$sampleID <- pca.names$V1

df.admix.k3 <- df.admix.k3 %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -sampleID) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.admix.k3$loc <- factor(ifelse(grepl("ACA", df.admix.k3$sampleID), "American", "European"))
df.admix.k3$pop <- as.character(gsub('[[:digit:]]+', '', df.admix.k3$sampleID))
pop.names <- read.csv("data/metadata/eel_pop_names.csv", stringsAsFactors = F)
df.admix.k3 <- df.admix.k3 %>% left_join(pop.names, by = "pop")
df.admix.k3$Location <- factor(df.admix.k3$Location, levels = location.order)

df.admix.k3$Location <- factor(df.admix.k3$Location, levels = location.order)
df.admix.k3 <- df.admix.k3 %>% filter(!sampleID %in% deviant.list)
#https://luisdva.github.io/rstats/model-cluster-plots/

p.ngsadmix.k3 <- ggplot(df.admix.k3, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  #facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("grey", "#ef8a62", "#67a9cf"), guide = F) +
  scale_color_manual(values = c("grey", "#ef8a62", "#67a9cf"), guide = F)

p.ngsadmix.k3
ggsave("output/ngsadmix/admixture_with_canada_whole_genome_k3.png", width = 10, height = 4)


p.ngsadmix.k2.loc / p.ngsadmix.k3 / p.minor.ancestry + plot_annotation(tag_levels = 'A')
ggsave("output/ngsadmix/compare_k2_k3_minor_ancestry.pdf", width = 11, height = 8.5)


# without canada ----------------------------------------------------------

# NGSadmix method ---------------------------------------------------------
#files say super 1 but are whole genome

S.admix <- read.table("data/ngsadmix/subset_whole_genome/ngs_admix_k2_no_canada.qopt") # Reads results from selection scan
pca.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 
pca.names.nc <- pca.names %>% filter(!grepl("ACA", V1))
deviant.list <- c("AMAR7", "ATU6", "APM45", "AMAR10")

df.admix.k2 <- as.data.frame(S.admix)
df.admix.k2$sampleID <- pca.names.nc$V1

df.admix.k2 <- df.admix.k2 %>% pivot_longer(names_to = "popGroup", values_to = "prob", cols = -sampleID) %>% 
  mutate(popGroup = gsub("V", "", popGroup))

df.admix.k2$loc <- factor(ifelse(grepl("ACA", df.admix.k2$sampleID), "American", "European"))
df.admix.k2$pop <- as.character(gsub('[[:digit:]]+', '', df.admix.k2$sampleID))
pop.names <- read.csv("data/metadata/eel_pop_names.csv", stringsAsFactors = F)
df.admix.k2 <- df.admix.k2 %>% left_join(pop.names, by = "pop")

df.admix.k2 <- df.admix.k2 %>% filter(!sampleID %in% deviant.list)
#https://luisdva.github.io/rstats/model-cluster-plots/
df.admix.k2$Location <- factor(df.admix.k2$Location, levels = location.order)


p.ngsadmix.no.canda <- ggplot(df.admix.k2, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(aes(color = popGroup), size = 0.1) +
  #facet_grid(~loc, switch = "x", scales = "free", space = "free") +
  facet_grid(~Location, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE) + xlab(NULL) +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf"), guide = F)+
  scale_color_manual(values = c("#ef8a62", "#67a9cf"), guide = F)
p.ngsadmix.no.canda
ggsave("output/ngsadmix/admixture_with_canada_whole_genome_k2_no_can.png", width = 10, height = 4)
ggsave("output/ngsadmix/admixture_with_canada_whole_genome_k2_no_can.pdf", width = 10, height = 4)


p.ngsadmix.k2.loc / p.ngsadmix.k3 / p.ngsadmix.no.canda/ p.minor.ancestry + plot_annotation(tag_levels = 'A')
ggsave("output/ngsadmix/compare_k2_k3_nc_minor_ancestry.pdf", width = 8.5, height = 11)
ggsave("output/ngsadmix/compare_k2_k3_nc_minor_ancestry.png", width = 8.5, height = 11)


