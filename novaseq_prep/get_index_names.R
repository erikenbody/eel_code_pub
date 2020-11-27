library(dplyr)
library(tidyr)
library(Biostrings)
# Pull in indexes and make some helpful files-------------------------------------------------------------------------

idx.layout <- read.csv("data/novaseq_data_input/Index_primer_layout.csv", header = T)
idx.layout <- idx.layout %>% 
  mutate(plate.idx = paste(Plate, Sample_Well, sep = "-"))

#I found that local run manager only accepts indexs with a prefix from the plate in the illumina kit. Found the identity of these here:
#http://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nextera-xt/nextera-xt-library-prep-reference-guide-15031942-05.pdf
idx.layout$I7_MiniSeq_ID <- ifelse(gsub("N","",idx.layout$I7_Index_ID) <= 715, paste(idx.layout$I7_Index_ID, "A", sep = "-"), paste(idx.layout$I7_Index_ID, "B",sep = "-"))
idx.layout$I5_MiniSeq_ID <- ifelse(gsub("S","",idx.layout$I5_Index_ID) <= 511, paste(idx.layout$I5_Index_ID, "A", sep = "-"), paste(idx.layout$I5_Index_ID, "C",sep = "-"))

i7 <- idx.layout %>% select(I7_Index_ID, I7_MiniSeq_ID, index) %>% group_by(index) %>% filter(!duplicated(index)) %>% mutate(IndexReadNumber = 1)
i5 <- idx.layout %>% select(I5_Index_ID, I5_MiniSeq_ID, index2) %>% group_by(index2) %>% filter(!duplicated(index2)) %>% mutate(IndexReadNumber = 2)
write.csv(i7, "output/i7.csv", row.names = F)
write.csv(i5, "output/i5.csv", row.names = F)


# make miniseq file -------------------------------------------------------
#removed greece
pools.idx <- read.csv("data/novaseq_data_input/Eel pilot project plate records.csv", header = T)
head(pools.idx)
#format sample wells to match index dataframe
pools.idx$Sample <- as.character(pools.idx$Sample)
pools.idx$Sample_Well <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", pools.idx$Pos, perl = TRUE)
pools.idx <- pools.idx %>% 
  mutate(plate.idx = paste(Index_Plate, Sample_Well, sep = "-")) %>% 
  select(Sample, plate.idx, Pool)

idx.pool <- left_join(pools.idx, idx.layout, by = "plate.idx") 

# make novaseq file -------------------------------------------------------

idx.pool$index2_4nova <- sapply(idx.pool$index2, function(x) as.character(reverseComplement(DNAString(x))))
novaseq.idx.pool <- idx.pool %>% select(Sample, index, index2_4nova, Pool) %>% 
  mutate(CUSTOM_INDEX = paste(index, index2_4nova, sep = "-")) %>% 
  select(Pool, Sample, CUSTOM_INDEX, -index, -index2_4nova) %>% 
  dplyr::rename(LIBRARY_ID = Sample)

write.csv(novaseq.idx.pool, "data/novaseq_data_input/eel_sample_indexes.csv", row.names = F, quote = F)
