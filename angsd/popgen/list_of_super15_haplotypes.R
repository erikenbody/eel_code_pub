library(tidyverse)
load("data/pcangsd/individual_group_assignment.rdata")
head(Beagle_ind_list)
table(Beagle_ind_list$super_15_clust_v2)

for (group in unique(Beagle_ind_list$super_15_clust_v2)){
  Beagle_ind_list %>% filter(super_15_clust_v2 == group) %>% 
    select(V1) %>% 
    write.table(paste0("data/pcangsd/super_15_",group,".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
}
