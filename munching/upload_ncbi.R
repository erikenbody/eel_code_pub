ncbi.names <- read.table("data/pcangsd/eel_names.txt", stringsAsFactors = F) #custom formatted 
pop.names <- read.csv("data/metadata/metadata_4_ncbi.csv", stringsAsFactors = F)
names(ncbi.names) <- "sample_name"
ncbi.names$pop <- as.character(gsub('[[:digit:]]+', '', ncbi.names$sample))
ncbi.names <- ncbi.names %>% left_join(pop.names, by = "pop")
ncbi.names$organism <- ncbi.names$Species
ncbi.names$isolate <- ncbi.names$pop
ncbi.names$dev_stage <- "young"
ncbi.names$sex <- "missing"
ncbi.names$tissue <- "not applicable"
ncbi.names$biomaterial_provider <- "Leif Andersson"
ncbi.names$collected_by <- "Hakan Wickstrom"
ncbi.names$collection_date <- ncbi.names$Date
ncbi.names$geo_loc_name <- paste(ncbi.names$Location, ncbi.names$Detail, sep=":")
ncbi.names$lat_lon <- paste(ncbi.names$Lattitude, ncbi.names$Longitude, sep = " ")
ncbi.names$sample_type <- "tissue sample"

ncbi.names %>% select(sample_name, organism, isolate, dev_stage, sex, tissue, biomaterial_provider, collected_by, collection_date, geo_loc_name, lat_lon, sample_type) %>% 
  write_csv("docs/SRA_SUB/eel_metadata.csv")


df.fastq <- read.table("docs/SRA_SUB/R1_names.txt")
df.fastq.out <- df.fastq %>% mutate(name1 = gsub("TD-2544-","", V1)) %>% 
  separate(name1, into =c("sample",NA,NA,NA,NA)) %>% 
  mutate(R2 = gsub("R1","R2", V1)) %>% 
  dplyr::rename(R1 = V1)

write.csv(df.fastq.out, "docs/SRA_SUB/fastq_names_out.csv")

           