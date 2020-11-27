library(tidyverse)
library(patchwork)
library(zoo)
library(data.table)
# -------------------------------------------------------------------------

l.theta <- list()
l.watt <- list()
l.sum <- list()
l.dfs <- list()
for (pop in file.list){
  df.theta <- read.table(pop)
  pop.name <- gsub(".theta.thetasWindow.gz.pestPG","",basename(pop))
  l.theta[[pop.name]] <- mean((df.theta$V5 / df.theta$V14), na.rm = T)
  l.watt[[pop.name]] <- mean((df.theta$V4 /df.theta$V14), na.rm = T)
  l.sum[[pop.name]] <- sum(df.theta$V14, na.rm = T)
  df.theta$pop <- pop.name
  l.dfs[[pop.name]] <- df.theta
}

df.x2 <- data.frame(pop = names(l.theta),
                   pairwise.nuc = round(unlist(l.theta[1:11]), 4),
                   wattersons.theta = round(unlist(l.watt[1:11]), 4),
                   number.sites = round(unlist(l.sum[1:11]), 4))

write.csv(df.x2, "output/thetas/pop_tP_mean_unfolded_nSites.csv", row.names = F)
mean.theta.nsites <- df.x2
df.theta.nsites <- bind_rows(l.dfs)
names(df.theta.nsites) <- c("window", "Chr", "WinCenter","tW", "tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites", "pop")


# cleaner summary ---------------------------------------------------------
#always filter windows > 10 sites

clean.summary.thetas <- df.theta.nsites %>% filter(nSites > 10) %>%
  group_by(pop) %>%
  summarise(mean.pi = round(mean(tP / nSites),3),
            sd.pi = round(sd(tP / nSites),2),
            mean.wat = round(mean(tW / nSites),3),
            sd.wat = round(sd(tW / nSites),2)) %>%
  mutate(out.pi = paste0(mean.pi,"±",sd.pi),
         out.wat = paste0(mean.wat,"±",sd.wat),)
write.csv(clean.summary.thetas, "output/thetas/clean_summaries/clean_theta_folded_perpop.csv")


clean.species.thetas <- df.theta.nsites %>% filter(nSites > 10) %>%
  mutate(species = ifelse(pop == "ACA", "American", "Euro")) %>%
  group_by(species) %>%
  summarise(mean.pi = round(mean(tP / nSites),3),
            sd.pi = round(sd(tP / nSites),2),
            mean.wat = round(mean(tW / nSites),3),
            sd.wat = round(sd(tW / nSites),2)) %>%
  mutate(out.pi = paste0(mean.pi, "±",sd.pi),
         out.wat = paste0(mean.wat, "±",sd.wat),)
write.csv(clean.species.thetas, "output/thetas/clean_summaries/clean_theta_folded_species.csv")
