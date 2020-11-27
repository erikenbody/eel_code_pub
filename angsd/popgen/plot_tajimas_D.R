library(tidyverse)
library(patchwork)

# creat ind list ----------------------------------------------------------

load("data/pcangsd/individual_group_assignment.rdata")
head(Beagle_ind_list)
table(Beagle_ind_list$super_15_clust_v2)

for (group in unique(Beagle_ind_list$super_15_clust_v2)){
  Beagle_ind_list %>% filter(super_15_clust_v2 == group) %>% 
    select(V1) %>% 
    write.table(paste0("data/pcangsd/super_15_",group,".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
}

#ran in angsd with only filters not MAF

# -------------------------------------------------------------------------

df.sel15 <- read.csv("output/pcangsd/selection_SUPER_15.csv")
df.sel15 <- df.sel15 %>% separate(pos, into = c(NA,"CHR","position"))
df.sel15 <- df.sel15 %>% mutate(position = as.numeric(position))

df.theta1 <- read.table("data/thetas/tajimas_d/group1_SUPER_15.theta.thetasWindow.gz.pestPG", comment.char = "", header = T)[,-1]
df.theta1$Tajima.ps <- df.theta1$Tajima / 10000
df.theta1$Tajima.rm <- zoo::rollmean(df.theta1$Tajima.ps,5,fill=NA)

df.theta4 <- read.table("data/thetas/tajimas_d/group4_SUPER_15.theta.thetasWindow.gz.pestPG", comment.char = "", header = T)[,-1]
df.theta4$Tajima.ps <- df.theta4$Tajima / 10000 
df.theta4$Tajima.rm <- zoo::rollmean(df.theta4$Tajima.ps,5,fill=NA)

p.sel <- ggplot()+
  geom_point(data = subset(df.sel15, position > 10000000 & position < 14000000), aes(x = position/1000000, y = -log(pc1.pval)), color = "black", alpha = 0.5, size = 1.5) +
  geom_line(data = subset(df.sel15, position > 10000000 & position < 14000000), aes(x = position/1000000, y = -log(pc1.pval.rollmean)), color = "red", size = 1.5) +
  theme_bw() + 
  theme(text=element_text(size=21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("-log(PC1:pvalue)") + xlab(NULL)

p.taj <- ggplot() + 
  geom_point(data = subset(df.theta1, WinCenter > 10000000 & WinCenter < 14000000), aes(x = WinCenter/1000000, y = (Tajima.ps*10^3)), size = 0.9, color = "red", alpha = 0.8) +
  geom_point(data = subset(df.theta4, WinCenter > 10000000 & WinCenter < 14000000), aes(x = WinCenter/1000000, y = (Tajima.ps*10^3)), size = 0.9, color = "blue", alpha = 0.8) +
  geom_line(data = subset(df.theta1, WinCenter >  10000000 & WinCenter < 14000000), aes(x = WinCenter/1000000,  y = (Tajima.rm*10^3)), color = "red", size = 1.5) +
  geom_line(data = subset(df.theta4, WinCenter >  10000000 & WinCenter < 14000000), aes(x = WinCenter/1000000,  y = (Tajima.rm*10^3)), color = "blue", size = 1.5) +
  theme_bw() +xlab("Mbp") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text=element_text(size=21)) +
  labs(y="Tajima's D * 10^3")

p.sel / p.taj
ggsave("output/thetas/tajimas_D_SUPER_15.png", height = 10, width = 8,dpi = 320)


# 1kb windows -------------------------------------------------------------

df.theta1 <- read.table("data/thetas/tajimas_d/group1_SUPER_15.theta.1kb.thetasWindow.gz.pestPG", comment.char = "", header = T)[,-1]
df.theta1$Tajima.ps <- df.theta1$Tajima / 1000
df.theta1$Tajima.rm <- zoo::rollmean(df.theta1$Tajima.ps,5,fill=NA)

df.theta4 <- read.table("data/thetas/tajimas_d/group4_SUPER_15.theta.1kb.thetasWindow.gz.pestPG", comment.char = "", header = T)[,-1]
df.theta4$Tajima.ps <- df.theta4$Tajima / 1000 
df.theta4$Tajima.rm <- zoo::rollmean(df.theta4$Tajima.ps,5,fill=NA)


p.sel2 <- ggplot()+
  geom_point(data = subset(df.sel15, position > 11900000 & position < 12200000), aes(x = position/1000000, y = -log(pc1.pval)), color = "black", alpha = 0.5, size = 1.5) +
  geom_line(data = subset(df.sel15, position > 11900000 & position < 12200000), aes(x = position/1000000, y = -log(pc1.pval.rollmean)), color = "red", size = 1.5) +
  theme_bw() + 
  theme(text=element_text(size=21)) +
  ylab("-log(PC1:pvalue)") + xlab("Mbp")

p.taj2 <- ggplot() + 
  geom_point(data = subset(df.theta1, WinCenter > 11900000 & WinCenter < 12200000), aes(x = WinCenter/1000000, y = (Tajima.ps*10^3)), size = 0.9, color = "red", alpha = 0.8) +
  geom_point(data = subset(df.theta4, WinCenter > 11900000 & WinCenter < 12200000), aes(x = WinCenter/1000000, y = (Tajima.ps*10^3)), size = 0.9, color = "blue", alpha = 0.8) +
  geom_line(data = subset(df.theta1, WinCenter >  11900000 & WinCenter < 12200000), aes(x = WinCenter/1000000,  y = (Tajima.rm*10^3)), color = "red", size = 1) +
  geom_line(data = subset(df.theta4, WinCenter >  11900000 & WinCenter < 12200000), aes(x = WinCenter/1000000,  y = (Tajima.rm*10^3)), color = "blue", size = 1) +
  theme_bw() +xlab("Mbp") + 
  theme(text=element_text(size=21)) +
  labs(y="Tajima's D * 10^3")

p.sel2 / p.taj2
ggsave("output/thetas/tajimas_D_SUPER_15_1kb.png", height = 10, width = 8)
