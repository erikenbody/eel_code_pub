#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Eel data processing
R1_files <- read.table("~/Projects/Eel/data/mappings/R1_files.txt", stringsAsFactors = F)
R1_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_001.fastq.gz", "\\1", R1_files[,1])
R1_files[,"Short_ID"] <- sub("_S[0-9]+_L[0-9]+", "", R1_files[,"ID"])
R2_files <- read.table("~/Projects/Eel/data/mappings/R2_files.txt", stringsAsFactors = F)
R2_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_001.fastq.gz", "\\1", R2_files[,1])
all(R2_files[,"ID"] == R1_files[,"ID"])
#Already in correct order, no sorting needed

require(Biostrings)
#Reference genome
fAngAng1_pri <- readDNAStringSet("~/Projects/Eel/data/reference/fAngAng1.pri.cur.20200204.fasta.gz")
write(x = names(fAngAng1_pri), file = "~/Projects/Eel/data/reference/fAngAng1_pri_chromosomes.txt")

#List of bam files
#find /proj/uppstore2017191/private/users/mats/Eel/bams/ -maxdepth 5 -type f -name '*.bam' > eel_bams.txt

#One bam seems to be missing
eel_bams <- read.table("~/Projects/Eel/data/mappings/eel_bams.txt", stringsAsFactors = F)
eel_bams[,"ID"] <- sub(".+/Sample_([A-Za-z0-9_-]+).sort.MarkDup.bam", "\\1", eel_bams[,1])
which(!(R1_files$Short_ID %in% eel_bams$ID))
#[1] 471 - Checked slurm-logs, ran out of time

#Complete version
#find /proj/uppstore2017191/private/users/mats/Eel/bams/ -maxdepth 5 -type f -name '*.bam' > eel_bams_v2.txt
#eel_bams <- read.table("~/Projects/Eel/data/mappings/eel_bams_v2.txt", stringsAsFactors = F)
#eel_bams[,"ID"] <- sub(".+/Sample_([A-Za-z0-9_-]+).sort.MarkDup.bam", "\\1", eel_bams[,1])
#eel_bams[,"group"] <- sub("TD-2544-([A-Za-z]+)[0-9]+", "\\1", eel_bams[,"ID"])

#Update after re-locating to new storage project
#find /proj/snic2020-2-19/private/eel/alignments -maxdepth 5 -type f -name '*.bam' > eel_bams_re_loc.txt
eel_bams <- read.table("~/Projects/Eel/data/bam_re_loc/eel_bams_re_loc.txt", stringsAsFactors = F)
eel_bams[,"ID"] <- sub(".+/Sample_([A-Za-z0-9_-]+).sort.MarkDup.bam", "\\1", eel_bams[,1])
eel_bams[,"group"] <- sub("TD-2544-([A-Za-z]+)[0-9]+", "\\1", eel_bams[,"ID"])


#Making links for SRA uploads
cat R1_files.txt | while read r1_file; do ln -s $r1_file; done
cat R2_files.txt | while read r2_file; do ln -s $r2_file; done

#[matsp@r187 reads]$ lftp subftp@ftp-private.ncbi.nlm.nih.gov
#Password: 
#  lftp subftp@ftp-private.ncbi.nlm.nih.gov:~> cd uploads/mats.pettersson_imbim.uu.se_7UxcCkRi
#cd ok, cwd=/uploads/mats.pettersson_imbim.uu.se_7UxcCkRi                     
#lftp subftp@ftp-private.ncbi.nlm.nih.gov:/uploads/mats.pettersson_imbim.uu.se_7UxcCkRi> mkdir eel_low_pass
#mkdir ok, `eel_low_pass' created                
#lftp subftp@ftp-private.ncbi.nlm.nih.gov:/uploads/mats.pettersson_imbim.uu.se_7UxcCkRi> cd eel_low_pass/
#lftp subftp@ftp-private.ncbi.nlm.nih.gov:/uploads/mats.pettersson_imbim.uu.se_7UxcCkRi/eel_low_pass> mput -c ./*.fastq.gz


#Trying some possible groupings
#grep -E "ASRI|SM|AST" eel_bams_v2.txt > swedish_bams.txt 
#grep -E "AAF|AI|AES" eel_bams_v2.txt > cont_europe_bams.txt 

#Per-group lists
for(g_name in unique(eel_bams[,"group"])){
  #write(file = paste("~/Projects/Eel/data/mappings/", g_name, "_bams.txt", sep = ""), x = eel_bams[eel_bams[,"group"] == g_name,1], ncolumns = 1)
  write(file = paste("~/Projects/Eel/data/bam_re_loc/", g_name, "_bams.txt", sep = ""), x = eel_bams[eel_bams[,"group"] == g_name,1], ncolumns = 1)
}

#Checking example chr - SUPER_19
SUPER_19_maf <- read.table("~/Projects/Eel/data/All_MnM_SUPER_19.mafs.gz", header = T, stringsAsFactors = F, sep = "\t")

#Angsd site file (CHR\tPOS)
#awk '{print $1"\t"$2}' All_MnM_SUPER_19.mafs > SUPER_19_maf_0.05_sites.txt
#/home/matsp/private/Software/angsd/angsd/angsd sites index SUPER_19_maf_0.05_sites.txt

SUPER_19_maf[,"SNP_ID"] <- paste(SUPER_19_maf$chromo,SUPER_19_maf$position, sep = "_")

SUPER_19_maf_SWE <- read.table("~/Projects/Eel/data/SWE_MnM_SUPER_19.mafs.gz", header = T, stringsAsFactors = F, sep = "\t")
SUPER_19_maf_SWE[,"SNP_ID"] <- paste(SUPER_19_maf_SWE$chromo,SUPER_19_maf_SWE$position, sep = "_")

SUPER_19_maf_EU <- read.table("~/Projects/Eel/data/EU_MnM_SUPER_19.mafs.gz", header = T, stringsAsFactors = F, sep = "\t")
SUPER_19_maf_EU[,"SNP_ID"] <- paste(SUPER_19_maf_EU$chromo,SUPER_19_maf_EU$position, sep = "_")

SUPER_19_maf[,c("SWE_maf", "SWE_ID")] <- SUPER_19_maf_SWE[match(SUPER_19_maf$SNP_ID, SUPER_19_maf_SWE$SNP_ID), c("knownEM", "SNP_ID")]
SUPER_19_maf[,c("EU_maf", "EU_ID")] <- SUPER_19_maf_EU[match(SUPER_19_maf$SNP_ID, SUPER_19_maf_EU$SNP_ID), c("knownEM", "SNP_ID")] 

SUPER_19_complete <- !is.na(SUPER_19_maf$SWE_maf) & !is.na(SUPER_19_maf$EU_maf)

SUPER_19_maf <- SUPER_19_maf[SUPER_19_complete,]
sum(SUPER_19_maf$SNP_ID == SUPER_19_maf$EU_ID & SUPER_19_maf$SNP_ID == SUPER_19_maf$SWE_ID) #Checks out
SUPER_19_maf[,"DAF"] <- abs(SUPER_19_maf$SWE_maf - SUPER_19_maf$EU_maf)

#Visualisation
ind_co <- 400
png("~/Projects/Eel/doc/SUPER_19_SWEvEU.png", width = 1000)
plot(x = SUPER_19_maf$position[SUPER_19_maf$nInd > ind_co], y = SUPER_19_maf$DAF[SUPER_19_maf$nInd > ind_co], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = "SUPER_19", col = "grey70")
rolling_DAF <- filter(x=SUPER_19_maf$DAF[SUPER_19_maf$nInd > ind_co], rep(1/100, 100))
points(x = SUPER_19_maf$position[SUPER_19_maf$nInd > ind_co], y = rolling_DAF, pch = 16, cex = 0.5, ylab = "DAF (100 SNP average)", xlab = "Position", main = "SUPER_19", col = "firebrick")
dev.off()

SUPER_19_maf <- read_and_plot_contrast("SUPER_19")
SUPER_1_maf <- read_and_plot_contrast("SUPER_1")

for(t_chr in paste0("SUPER_", c(2:9, 11:18))){
  read_and_plot_contrast(t_chr)
}

#Regions of interest
read_and_plot_contrast("SUPER_3", x_lim = c(3.2e7,3.25e7))
read_and_plot_contrast("SUPER_3", x_lim = c(3.22e7,3.24e7))
read_and_plot_contrast("SUPER_1", x_lim = c(8.1e7,8.3e7))
read_and_plot_contrast("SUPER_1", x_lim = c(8.11e7,8.13e7))
read_and_plot_contrast("SUPER_12", x_lim = c(1.95e7,2.05e7))
read_and_plot_contrast("SUPER_12", x_lim = c(1.98e7,2.00e7))
read_and_plot_contrast("SUPER_16", x_lim = c(2.00e7,2.24e7))
read_and_plot_contrast("SUPER_16", x_lim = c(2.05e7,2.15e7))
read_and_plot_contrast("SUPER_16", x_lim = c(2.15e7,2.25e7))
read_and_plot_contrast("SUPER_16", x_lim = c(2.16e7,2.19e7))


for(t_super in paste0("SUPER_", c(1:19))){
  assign(paste0(t_super,"_freq"), read_maf_files(t_super, g_ids))
  save(list = paste0(t_super,"_freq"), file = paste0("~/Projects/Eel/data/all_maf/", t_super, "_freq.RData"))
}

#5% cut-off version
for(t_super in paste0("SUPER_", c(1:17))){
  assign(paste0(t_super,"_0.05_freq"), read_maf_files(t_super, g_ids, total_dir = "~/Projects/Eel/data/all_maf/maf_co_0.05/", per_group_dir = "~/Projects/Eel/data/per_group_maf/maf_co_0.05/"))
  save(list = paste0(t_super,"_0.05_freq"), file = paste0("~/Projects/Eel/data/all_maf/maf_co_0.05/", t_super, "_0.05_freq.RData"))
}

#Test example
SUPER_16_22MB_freq <- SUPER_16_freq[SUPER_16_freq$position > 2.16e7 & SUPER_16_freq$position < 2.19e7 & SUPER_16_freq$nInd > 300,-(1:8)]
rownames(SUPER_16_22MB_freq) <- SUPER_16_freq[SUPER_16_freq$position > 2.16e7 & SUPER_16_freq$position < 2.19e7 & SUPER_16_freq$nInd > 300,"position"]
high_daf <- abs(rowMeans(SUPER_16_22MB_freq[,c("AAF", "AI", "AES")]) - rowMeans(SUPER_16_22MB_freq[,c("ASRI", "SMS", "AST")])) > 0.2
heatmap(t(as.matrix(SUPER_16_22MB_freq[high_daf,])), scale = "none", Colv = NA)
plot(y = abs(rowMeans(SUPER_16_22MB_freq[,c("AAF", "AES")]) - rowMeans(SUPER_16_22MB_freq[,c("ASRI", "SMS", "AST")])), x = rownames(SUPER_16_22MB_freq), pch = 16, xlab = "Position on SUPER_16", main = "AES + AAF v SWE", ylab = "DAF")

#Regions of interest
g_ids <- unique(eel_bams[,"group"])

#SUPER_12_freq <- read_maf_files("SUPER_12", g_ids )
#save(SUPER_12_freq, file = "~/Projects/Eel/data/all_maf/SUPER_12_freq.RData")

SUPER_12_20MB_GR <- GRanges(seqnames = "SUPER_12", ranges = IRanges(start =  1.98e7, end = 2.00e7))
plot_region_hm(SUPER_12_20MB_GR, g1_ids = midAtl_ids, g2_ids = non_midAtl_ids, daf_co = 0.08)


#SUPER_1_freq <- read_maf_files("SUPER_1", g_ids)
#save(SUPER_1_freq, file = "~/Projects/Eel/data/all_maf/SUPER_1_freq.RData")
SUPER_1_81MB_GR <- GRanges(seqnames = "SUPER_1", ranges = IRanges(start = 8.118e7, end = 8.120e7))
#plot_region_hm(SUPER_1_81MB_GR, g1_ids = midAtl_ids, g2_ids = non_midAtl_ids, daf_co = 0.12)
plot_region_hm(SUPER_1_81MB_GR, g1_ids = midAtl_ids, g2_ids = non_midAtl_ids, daf_co = 0.12, plot_dir = "~/Projects/Eel/doc/MidAtl_0.05/", target_freq = SUPER_1_0.05_freq)

SUPER_3_10MB_GR <- GRanges(seqnames = "SUPER_3", ranges = IRanges(start =  1.02e7, end = 1.025e7))
plot_region_hm(SUPER_3_10MB_GR, g1_ids = midAtl_ids, g2_ids = non_midAtl_ids)

SUPER_16_15MB_GR <- GRanges(seqnames = "SUPER_16", ranges = IRanges(start =  1.453e7, end = 1.46e7))
plot_region_hm(SUPER_16_15MB_GR, g1_ids = midAtl_ids, g2_ids = non_midAtl_ids)



baltic_ids <- c("SMS","AST", "ALI")
non_baltic_ids <- g_ids[!g_ids %in% c(baltic_ids, "ACA")]
#plot_contrast(g1_ids = baltic_ids, g2_ids = non_baltic_ids, contrast_folder = "baltic_v_rest", chr_list = 1:19)
baltic_DAF <- plot_contrast(g1_ids = baltic_ids, g2_ids = non_baltic_ids, contrast_folder = "~/Projects/Eel/doc/Baltic_0.05_rerun/", chr_list = 1:19)
summary(baltic_DAF$DAF)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.009511 0.021159 0.027122 0.038259 0.296586 

midAtl_ids <- c("AAF","AES", "AI")
non_midAtl_ids <- g_ids[!g_ids %in% c(midAtl_ids, "ACA")]
#plot_contrast(g1_ids = midAtl_ids, g2_ids = non_midAtl_ids, contrast_folder = "MidAtl_v_rest", chr_list = 1:19)
#midAtl_DAF <- 
midAtl_DAF <- plot_contrast(g1_ids = midAtl_ids, g2_ids = non_midAtl_ids, contrast_folder = "~/Projects/Eel/doc/MidAtl_0.05_rerun_v2/", chr_list = 1:19)
summary(midAtl_DAF$DAF)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.009656 0.021164 0.027028 0.038093 0.296068 

southAtl_ids <- c("AMAR","APM")
non_southAtl_ids <- g_ids[!g_ids %in% c(southAtl_ids, "ACA")]
plot_contrast(g1_ids = southAtl_ids, g2_ids = non_southAtl_ids, contrast_folder = "SouthAtl_v_rest", chr_list = 1:19)

medSea_ids <- c("ATU")
non_medSea_ids  <- g_ids[!g_ids %in% c(medSea_ids, "ACA")]
plot_contrast(g1_ids = medSea_ids, g2_ids = non_medSea_ids, contrast_folder = "MedSea_v_rest", chr_list = 16:19)

american_ids <- c("ACA")
non_american_ids  <- g_ids[!g_ids %in% c("ACA")]
#plot_contrast(g1_ids = american_ids, g2_ids = non_american_ids, contrast_folder = "Canada_v_rest", chr_list = 1:19)
ACA_daf <- plot_contrast(g1_ids = american_ids, g2_ids = non_american_ids, contrast_folder = "~/Projects/Eel/doc/ACA_rerun/", chr_list = 1:19)


motala_ids <- c("SMS")
non_motala_ids  <- g_ids[!g_ids %in% c("ACA", "SMS")]
plot_contrast(g1_ids = motala_ids, g2_ids = non_motala_ids, contrast_folder = "~/Projects/Eel/doc/Motala_0.05/", chr_list = 1:19)



#Individual genoytpe probabilies
#SUPER_15_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_19_SUPER_15.beagle.gz", header = T, sep = "\t", stringsAsFactors = F)
#Very resource heavy, using zgrep & awk to subset first

#Order of individuals in the beagle files
Beagle_ind_list <- read.table(file = "~/Projects/Eel/data/geno_likelihood/eel_bams_no_canada.txt", stringsAsFactors = F, header = F, sep = "\t")
Beagle_ind_list[,"ID"] <- sub(".+/Sample_TD-2544-([A-Za-z0-9_-]+).sort.MarkDup.bam", "\\1", Beagle_ind_list[,1])
Beagle_ind_list[,"GROUP"] <- sub("[0-9]+", "", Beagle_ind_list[,"ID"])


#Subset based on known region
SUPER_15_site_sel <- read.csv(file = "~/Projects/Eel/data/geno_likelihood/selection_SUPER_15.csv", stringsAsFactors = F)
SUPER_15_site_sel[,"position"] <- as.numeric(sub("SUPER_15_", "", SUPER_15_site_sel[,"pos"]))
#SUPER 15  peak 1
#SUPER_15_11870507 - SUPER_15_12246140
zgrep -m1 -n -o SUPER_15_11870507 SUPER_19_SUPER_15.beagle.gz
325503:SUPER_15_11870507
332246:SUPER_15_12246140
zcat SUPER_19_SUPER_15.beagle.gz | awk 'NR == 1 || (NR >=325503 && NR <= 332246) {print}' > SUPER_15_peak1.beagle 

SUPER_15_peak1_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak1.beagle.gz", header = T, sep = "\t", stringsAsFactors = F)
SUPER_15_peak1_GL[,"pos"] <- as.numeric(sub("SUPER_15_", "", SUPER_15_peak1_GL[,"marker"]))
SUPER_15_peak1_GL[,"pc1.pval"] <- SUPER_15_site_sel[match(SUPER_15_peak1_GL[,"pos"], SUPER_15_site_sel[,"position"]),"pc1.pval"]
SUPER_15_peak1_slide <- sliding_deviation(SUPER_15_peak1_GL, Beagle_ind_list$pc1_clust, pdf_file = "~/Projects/Eel/doc/SUPER_15_p1_slide.pdf", ref_ind = 1)


#SUPER 15 peak 2
#SUPER_15_33455012 - SUPER_15_33674630

#SUPER_13, start
SUPER_13_site_sel <- read.csv(file = "~/Projects/Eel/data/geno_likelihood/selection_SUPER_13.csv", stringsAsFactors = F)
SUPER_13_site_sel[,"position"] <- as.numeric(sub("SUPER_13_", "", SUPER_13_site_sel[,"pos"]))
which.min(abs(SUPER_13_site_sel$position - 1e7))# Approximate end of the sweep
SUPER_13_site_sel[360519,]# SUPER_13_9999957
zgrep -m1 -n -o SUPER_13_9999957 SUPER_19_SUPER_13.beagle.gz
363220:SUPER_13_9999957
zcat SUPER_19_SUPER_13.beagle.gz | awk 'NR <= 363220 {print}' > SUPER_13_block1.beagle
#Still very large, prnting names of informative markers
write(SUPER_13_site_sel$pos[SUPER_13_site_sel$position < 1e7 & SUPER_13_site_sel$pc1.pval < 1e-5], file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_block1_diagnostic.txt")
join <(sort SUPER_13_block1.beagle) <(sort SUPER_13_block1_diagnostic.txt) > SUPER_13_block1_diagnostic.beagle
#Patterns are odd, adding 10k random markers to get background behaviour
SUPER_13_diag <- SUPER_13_site_sel$pos[SUPER_13_site_sel$position < 1e7 & SUPER_13_site_sel$pc1.pval < 1e-5]
SUPER_13_bg <- SUPER_13_site_sel$pos[sample(which(SUPER_13_site_sel$position < 1e7), size = 1e4)]
write(unique(c(SUPER_13_diag, SUPER_13_bg)), file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_block1_ext_sample.txt")
join <(sort SUPER_13_block1.beagle) <(sort SUPER_13_block1_ext_sample.txt) > SUPER_13_ext_sample.beagle

SUPER_13_block1_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_block1_diagnostic.beagle", header = F, sep = " ", stringsAsFactors = F)
names(SUPER_13_block1_GL)[1:1338] <- names(SUPER_15_peak1_GL)[1:1338] #Failed to keep the header while subsetting
SUPER_13_block1_GL[,"pos"] <- as.numeric(sub("SUPER_13_", "", SUPER_13_block1_GL[,"marker"]))
SUPER_13_block1_GL <- SUPER_13_block1_GL[order(SUPER_13_block1_GL[,"pos"]),]
SUPER_13_block1_GL[,"pc1.pval"] <- SUPER_13_site_sel[match(SUPER_13_block1_GL[,"pos"], SUPER_13_site_sel[,"position"]),"pc1.pval"] 
SUPER_13_block1_kin <- region_kinship(SUPER_13_block1_GL, Beagle_ind_list)
save(SUPER_13_block1_kin, file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_block1_kinship.RData")
SUPER_13_block1_hclust <- hclust(as.dist(SUPER_13_block1_kin[[1]]))
pdf(file = "~/Projects/Eel/doc/SUPER_13_hm.pdf", width = 10, height = 10)
heatmap(SUPER_13_block1_kin[[1]], scale = "none", Rowv = as.dendrogram(SUPER_13_block1_hclust), as.dendrogram(SUPER_13_block1_hclust))
dev.off()



SUPER_13_block1_clust <- rect.hclust(SUPER_13_block1_hclust, k = 3)
for(i in 1:3) Beagle_ind_list[SUPER_13_block1_clust[[i]],"super_13_clust"] <- i
SUPER_13_block1_slide <- sliding_deviation(SUPER_13_block1_GL, Beagle_ind_list$super_13_clust, pdf_file = "~/Projects/Eel/doc/SUPER_13_slide.pdf", ref_ind = 41)
#ind 1 is group2, ind 2 is group 1 and ind 41 is group 3
plot(x = SUPER_13_block1_slide$freq[[1]], y = SUPER_13_block1_slide$freq[[2]], xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red")
abline(lm(SUPER_13_block1_slide$freq[[2]] ~ SUPER_13_block1_slide$freq[[1]]))
plot(x = SUPER_13_block1_slide$freq[[1]], y = SUPER_13_block1_slide$freq[[3]], xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red")

plot(x = SUPER_13_block1_GL$pos, y = abs(SUPER_13_block1_slide$freq[[3]]-SUPER_13_block1_slide$freq[[1]]), ylim = c(0,1))


#Patterns are wierd, trying random cluster assignment
rand_clust <- sample(Beagle_ind_list[,c("super_13_clust")])
SUPER_13_random_slide <- sliding_deviation(SUPER_13_block1_GL, rand_clust, pdf_file = "~/Projects/Eel/doc/SUPER_13_random.pdf", ref_ind = 41)

plot(x = SUPER_13_random_slide$freq[[1]], y = SUPER_13_random_slide$freq[[2]], xlim = c(0,1), ylim = c(0,1))
abline(lm(SUPER_13_random_slide$freq[[2]] ~ SUPER_13_random_slide$freq[[1]]))
abline(a = 0, b = 1, col = "red")

plot(x = SUPER_13_random_slide$freq[[1]], y = SUPER_13_random_slide$freq[[3]], xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red")
#Seems to check out, groups are very similar
table(Beagle_ind_list[,c("super_13_clust", "GROUP")]) #Groups 1 & 2 have strong association with some populations

#Repeating with expanded marker set
SUPER_13_ext_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_ext_sample.beagle", header = F, sep = " ", stringsAsFactors = F)
names(SUPER_13_ext_GL)[1:1338] <- names(SUPER_15_peak1_GL)[1:1338] #Failed to keep the header while subsetting
SUPER_13_ext_GL[,"pos"] <- as.numeric(sub("SUPER_13_", "", SUPER_13_ext_GL[,"marker"]))
SUPER_13_ext_GL <- SUPER_13_ext_GL[order(SUPER_13_ext_GL[,"pos"]),]
SUPER_13_ext_GL[,"pc1.pval"] <- SUPER_13_site_sel[match(SUPER_13_ext_GL[,"pos"], SUPER_13_site_sel[,"position"]),"pc1.pval"] 
diag_filter <- SUPER_13_ext_GL$marker %in% SUPER_13_diag
#Using clusters form the "diagnostic" set still
#SUPER_13_ext_slide <- sliding_deviation(SUPER_13_ext_GL, Beagle_ind_list$super_13_clust, pdf_file = "~/Projects/Eel/doc/SUPER_13_ext_slide.pdf", ref_ind = 41)

#Adding a pre-treatment step to set such sites to NA
SUPER_13_ext_GL_clean <- SUPER_13_ext_GL
MajorHomCols <- grep("Ind[0-9]+$", names(SUPER_13_ext_GL))
require(matrixStats)
for(i in 1:length(MajorHomCols)){
  ind_SDs <- rowSds(as.matrix(SUPER_13_ext_GL_clean[, MajorHomCols[i] + c(0,1,2)]), na.rm = T)
  SUPER_13_ext_GL_clean[ind_SDs == 0,  MajorHomCols[i] + c(0,1,2)] <- NA
}



SUPER_13_ext_kinship <- region_kinship(SUPER_13_ext_GL_clean[diag_filter,], Beagle_ind_list)
SUPER_13_ext_kinship[[1]][is.na(SUPER_13_ext_kinship[[1]])] <- mean(SUPER_13_ext_kinship[[1]], na.rm = T) #One pair APM45/AMAR10 has no value, setting it to the overall avearge
SUPER_13_ext_hclust <- hclust(as.dist(SUPER_13_ext_kinship[[1]]))
plot(SUPER_13_ext_hclust)
SUPER_13_ext_c_assign <- rect.hclust(SUPER_13_ext_hclust, k = 2)
for(i in 1:2) Beagle_ind_list[SUPER_13_ext_c_assign[[i]],"super_13_clust_v2"] <- i

SUPER_13_ext_freq <- kin_and_freq(SUPER_13_ext_GL_clean,Beagle_ind_list$super_13_clust_v2, ref_ind = 1) #Adjusted for NAs in the GL data and using refined clustering
save(SUPER_13_ext_freq, SUPER_13_ext_GL_clean, file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_ext.RData")

#pdf(file = "~/Projects/Eel/doc/SUPER_13_post_clean.pdf", width = 10, height = 10)
pdf(file = "~/Projects/Eel/doc/SUPER_13_post_clean_redraw.pdf", width = 10, height = 10)
#heatmap(SUPER_13_ext_kinship[[1]], scale = "none", Rowv = as.dendrogram(SUPER_13_ext_hclust), as.dendrogram(SUPER_13_ext_hclust))

plot(x = SUPER_13_ext_freq$type1_freq[!diag_filter], y = SUPER_13_ext_freq$type2_freq[!diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_13_ext_freq$type1_freq[diag_filter], y = SUPER_13_ext_freq$type2_freq[diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "darkorchid", cex = 0.3)

#Minor allele frequency versions
loci_to_flip <- SUPER_13_ext_freq$type1_freq >= 0.5
x_vec <- SUPER_13_ext_freq$type1_freq
x_vec[loci_to_flip] <- 1 - x_vec[loci_to_flip]
y_vec <- SUPER_13_ext_freq$type2_freq
y_vec[loci_to_flip] <- 1 - y_vec[loci_to_flip]
plot(x = x_vec[!diag_filter], y = y_vec[!diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3, xlab = "Frequecy in Type 1", ylab = "Frequecy in Type 2" )
points(x = x_vec[diag_filter], y = y_vec[diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "darkorchid", cex = 0.3)
points(x = x_vec[diag_filter], y = y_vec[diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "darkorchid", cex = 0.3)
abline(coef = c(0,1))

plot(x = SUPER_13_ext_GL_clean$pos[!diag_filter], y = abs(SUPER_13_ext_freq$type2_freq[!diag_filter]-SUPER_13_ext_freq$type1_freq[!diag_filter]), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_13_ext_GL_clean$pos[diag_filter], y = abs(SUPER_13_ext_freq$type2_freq[diag_filter]-SUPER_13_ext_freq$type1_freq[diag_filter]), pch = 16, col = "darkorchid", cex = 0.3)

dev.off()




#plot(x = SUPER_13_ext_slide$freq[[1]][!diag_filter], y = SUPER_13_ext_slide$freq[[2]][!diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
#points(x = SUPER_13_ext_slide$freq[[1]][diag_filter], y = SUPER_13_ext_slide$freq[[2]][diag_filter], col = "darkorchid", pch = 16, cex = 0.3)
#abline(lm(SUPER_13_ext_slide$freq[[2]][!diag_filter] ~ SUPER_13_ext_slide$freq[[1]][!diag_filter]))
#abline(lm(SUPER_13_ext_slide$freq[[2]][diag_filter] ~ SUPER_13_ext_slide$freq[[1]][diag_filter]), col = "darkorchid")
#abline(a = 0, b = 1, col = "red")

#plot(x = SUPER_13_ext_slide$freq[[1]][!diag_filter], y = SUPER_13_ext_slide$freq[[3]][!diag_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
#points(x = SUPER_13_ext_slide$freq[[1]][diag_filter], y = SUPER_13_ext_slide$freq[[3]][diag_filter], col = "darkorchid", pch = 16, cex = 0.3)
#abline(lm(SUPER_13_ext_slide$freq[[3]][!diag_filter] ~ SUPER_13_ext_slide$freq[[1]][!diag_filter]))
#abline(lm(SUPER_13_ext_slide$freq[[3]][diag_filter] ~ SUPER_13_ext_slide$freq[[1]][diag_filter]), col = "darkorchid")
#abline(a = 0, b = 1, col = "red")

#plot(x = SUPER_13_ext_GL[!diag_filter,"pos"] , y = abs(SUPER_13_ext_slide$freq[[3]]-SUPER_13_ext_slide$freq[[1]])[!diag_filter], ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
#points(x = SUPER_13_ext_GL[diag_filter,"pos"], y = abs(SUPER_13_ext_slide$freq[[3]]-SUPER_13_ext_slide$freq[[1]])[diag_filter], col = "darkorchid", pch = 16, cex = 0.3)

#plot(x = SUPER_13_ext_GL[!diag_filter,"pos"] , y = abs(SUPER_13_ext_slide$freq[[2]]-SUPER_13_ext_slide$freq[[1]])[!diag_filter], ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
#points(x = SUPER_13_ext_GL[diag_filter,"pos"], y = abs(SUPER_13_ext_slide$freq[[2]]-SUPER_13_ext_slide$freq[[1]])[diag_filter], col = "darkorchid", pch = 16, cex = 0.3)

#Checking marker information content in different subsamples
total_SDs <- rowSds(as.matrix(SUPER_13_ext_GL[, -(1:3)]))
g3_cols <- as.integer(outer(MajorHomCols[Beagle_ind_list$super_13_clust == 3], c(0,1,2), "+"))
g3_cols <- g3_cols[order(g3_cols)]
g2_cols <- as.integer(outer(MajorHomCols[Beagle_ind_list$super_13_clust == 2], c(0,1,2), "+"))
g2_cols <- g2_cols[order(g2_cols)]
g1_cols <- as.integer(outer(MajorHomCols[Beagle_ind_list$super_13_clust == 1], c(0,1,2), "+"))
g1_cols <- g1_cols[order(g1_cols)]

#Group-wise SD
g3_SDs <- rowSds(as.matrix(SUPER_13_ext_GL[, g3_cols]))

#Histograms
hist(unlist(SUPER_13_ext_GL[!diag_filter, g1_cols]))
hist(unlist(SUPER_13_ext_GL[!diag_filter, g1_cols]))
#Histograms indicate that this may be due to high number of unresolved sites (i.e. 1/3, 1/3, 1/3) in group 2





#SUPER_15_peak1, as an example
#Adding a pre-treatment step to set such sites to NA
SUPER_15_GL_clean <- SUPER_15_peak1_GL
MajorHomCols <- grep("Ind[0-9]+$", names(SUPER_15_peak1_GL))
pc1_filter <- SUPER_15_peak1_GL$pc1.pval < 1e-10
pc1_filter[is.na(pc1_filter)] <- F
for(i in 1:length(MajorHomCols)){
  ind_SDs <- rowSds(as.matrix(SUPER_15_GL_clean[, MajorHomCols[i] + c(0,1,2)]), na.rm = T)
  SUPER_15_GL_clean[ind_SDs == 0,  MajorHomCols[i] + c(0,1,2)] <- NA
}

#SUPER_15_freq <- kin_and_freq(SUPER_15_GL_clean, Beagle_ind_list$pc1_clust)


SUPER_15_kinship <- region_kinship(SUPER_15_GL_clean[pc1_filter,], Beagle_ind_list)
SUPER_15_kinship[[1]][is.na(SUPER_15_kinship[[1]])] <- mean(SUPER_15_kinship[[1]], na.rm = T) #Some pairs APM45/AMAR10 have no value, setting them to the overall avearge
SUPER_15_hclust <- hclust(as.dist(SUPER_15_kinship[[1]]))
plot(SUPER_15_hclust)
SUPER_15_c_assign <- rect.hclust(SUPER_15_hclust, k = 4)
for(i in 1:4) Beagle_ind_list[SUPER_15_c_assign[[i]],"super_15_clust_v2"] <- i
table(Beagle_ind_list[,c("pc1_clust", "super_15_clust_v2")])
heatmap(SUPER_15_kinship[[1]], scale = "none", Rowv = as.dendrogram(SUPER_15_hclust), as.dendrogram(SUPER_15_hclust))

SUPER_15_freq <- kin_and_freq(SUPER_15_GL_clean, Beagle_ind_list$super_15_clust_v2)
plot(x = SUPER_15_freq$type4_freq[!pc1_filter], y = SUPER_15_freq$type1_freq[!pc1_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_15_freq$type4_freq[pc1_filter], y = SUPER_15_freq$type1_freq[pc1_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "darkorchid", cex = 0.3)

plot(x = SUPER_15_GL_clean$pos[!pc1_filter], y = abs(SUPER_15_freq$type4_freq[!pc1_filter]-SUPER_15_freq$type1_freq[!pc1_filter]), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_15_GL_clean$pos[pc1_filter], y = abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]), pch = 16, col = "darkorchid", cex = 0.3)

#Output high DAF markers to intersect with snpEFF file
write(SUPER_15_GL_clean$pos[abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]) > 0.8], file = "~/Projects/Eel/data/annotation/SUPER_15_high_DAF.txt", ncolumns = 1)


#Below this line does NOT used data where uniformative markers are masked by NAs, view with caution but things seem mostly consistent

MajorHomCols <- grep("Ind[0-9]+$", names(SUPER_15_peak1_GL))
pc1_filter <- SUPER_15_peak1_GL$pc1.pval < 1e-10
pc1_filter[is.na(pc1_filter)] <- F

#kinship_mat <- matrix(nrow = length(Beagle_ind_list[,"ID"]), ncol = length(Beagle_ind_list[,"ID"]))
#rownames(kinship_mat) <- Beagle_ind_list[,"ID"]
#colnames(kinship_mat) <- Beagle_ind_list[,"ID"]
kinship_mat_high_pc1 <- matrix(nrow = length(Beagle_ind_list[,"ID"]), ncol = length(Beagle_ind_list[,"ID"]))
rownames(kinship_mat_high_pc1) <- Beagle_ind_list[,"ID"]
colnames(kinship_mat_high_pc1) <- Beagle_ind_list[,"ID"]
for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
  for(j in i:length(MajorHomCols)){ #length(MajorHomCols)
    #kinship_mat[i,j] <- sum(abs((2*SUPER_15_peak1_GL[,MajorHomCols[i]] + SUPER_15_peak1_GL[,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[,MajorHomCols[j]] + SUPER_15_peak1_GL[,MajorHomCols[j]+1])))/length(SUPER_15_peak1_GL[,MajorHomCols[i]])
    #kinship_mat[j,i] <- kinship_mat[i,j]
    kinship_mat_high_pc1[i,j] <- sum(abs((2*SUPER_15_peak1_GL[,MajorHomCols[i]] + SUPER_15_peak1_GL[,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[,MajorHomCols[j]] + SUPER_15_peak1_GL[,MajorHomCols[j]+1]))[pc1_filter])/sum(pc1_filter)
    kinship_mat_high_pc1[j,i] <- kinship_mat_high_pc1[i,j]
  }
}
save(kinship_mat, kinship_mat_high_pc1, file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak1_kinship.RData")
       
SUPER_15_peak1_hclust <- hclust(as.dist(kinship_mat))
#Using the the most infomative markers to divide teh indivduals
SUPER_15_peak1_pc1_hclust <- hclust(as.dist(kinship_mat_high_pc1))
plot(SUPER_15_peak1_pc1_hclust)
SUPER_15_peak1_pc1_rect <- rect.hclust(SUPER_15_peak1_pc1_hclust, k = 3)
for(i in 1:3) Beagle_ind_list[SUPER_15_peak1_pc1_rect[[i]],"pc1_clust"] <- i
 
write.table(table(Beagle_ind_list[,c("pc1_clust", "GROUP")]), file = "~/Projects/Eel/doc/SUPER_15_class_by_pop.txt", quote = F, sep = "\t")
#Average deviation for the classes
# kinship_avg_list <- list(t1_hom = numeric(sum(pc1_filter)), t2_hom = numeric(sum(pc1_filter)), het = numeric(sum(pc1_filter)))
# pc1_high_freq_avg_list <- list(t1_hom = numeric(sum(pc1_filter)), t2_hom = numeric(sum(pc1_filter)), het = numeric(sum(pc1_filter)))
# all_freq_avg_list <- list(t1_hom = numeric(dim(SUPER_15_peak1_GL)[1]), t2_hom = numeric(dim(SUPER_15_peak1_GL)[1]), het = numeric(dim(SUPER_15_peak1_GL)[1]))
# for(i in 2:length(MajorHomCols)){ #length(MajorHomCols)
#   kinship_vec <- abs((2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]]+ SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]+1]))
#   kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] <- kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + kinship_vec
#   pc1_high_freq_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] <- pc1_high_freq_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1])
#   all_freq_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] <- all_freq_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + (2*SUPER_15_peak1_GL[,MajorHomCols[i]] + SUPER_15_peak1_GL[,MajorHomCols[i]+1])
# }
# for(i in 1:3){
#   kinship_avg_list[[i]] <- kinship_avg_list[[i]]/sum(Beagle_ind_list[,"pc1_clust"] == i, na.rm = T)
#   pc1_high_freq_avg_list[[i]] <- pc1_high_freq_avg_list[[i]]/(2*sum(Beagle_ind_list[,"pc1_clust"] == i, na.rm = T))
#   all_freq_avg_list[[i]] <- all_freq_avg_list[[i]]/(2*sum(Beagle_ind_list[,"pc1_clust"] == i, na.rm = T))
# }

# pdf("~/Projects/Eel/doc/SUPER_15_peak1_kinship.pdf", width = 30)
# 
# plot(SUPER_15_peak1_hclust, main = "All SNPs", cex = 0.4)
# heatmap(kinship_mat, scale = "none", Rowv = as.dendrogram(SUPER_15_peak1_hclust), as.dendrogram(SUPER_15_peak1_hclust))
# 
# 
# plot(SUPER_15_peak1_pc1_hclust, main = "High pc1 SNPs", cex = 0.4)
# rect.hclust(SUPER_15_peak1_pc1_hclust, k = 3, border = c("blue", "red",  "darkorchid"))
# heatmap(kinship_mat_high_pc1, scale = "none", Rowv = as.dendrogram(SUPER_15_peak1_pc1_hclust), as.dendrogram(SUPER_15_peak1_pc1_hclust))
# 
# col_vec <- col2rgb(c("blue", "red",  "darkorchid"),alpha = T)
# col_vec["alpha",] <- 80
# 
# #Visualisation
# plot(x= SUPER_15_peak1_GL[pc1_filter,"pos"], y = SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]], type = "n", ylim = c(0,1.5), ylab = "Deviaiton from ind AAF10 (100 SNP rolling avg)", main = "SUPER_15", xlab = "Position")
# for(i in 2:length(MajorHomCols)){
#   kinship_vec <- abs((2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]]+ SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]+1]))
#   
#   lines(x = SUPER_15_peak1_GL[pc1_filter,"pos"], y = filter(kinship_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.1, col = rgb(t(col_vec[,Beagle_ind_list[i,"pc1_clust"]]), maxColorValue = 255)) 
# }
# plot(x= SUPER_15_peak1_GL[pc1_filter,"pos"], y = SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "SUPER_15", xlab = "Position")
# for(i in 2:length(MajorHomCols)){ #length(MajorHomCols)
#   kinship_dev_vec <- abs((2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]]+ SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]+1]))
#   kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + mean(kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]], na.rm = T)
#   #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
#   lines(x = SUPER_15_peak1_GL[pc1_filter,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.3, col = rgb(t(col_vec[,Beagle_ind_list[i,"pc1_clust"]]), maxColorValue = 255)) 
# }
# col_vec_hl <- col_vec
# col_vec_hl["alpha",] <- 150
# 
# for(j in 1:3){
#   plot(x= SUPER_15_peak1_GL[pc1_filter,"pos"], y = SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "SUPER_15", xlab = "Position")
#   for(i in which(Beagle_ind_list[,"pc1_clust"] != j)){ #length(MajorHomCols)
#     kinship_dev_vec <- abs((2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]]+ SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]+1]))
#     kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + mean(kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]], na.rm = T)
#     #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
#     if(i != 1) lines(x = SUPER_15_peak1_GL[pc1_filter,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.2, col = "grey70")
#     #Individual 1 is the comparison reference.
#   }
#   for(i in which(Beagle_ind_list[,"pc1_clust"] == j)){ #length(MajorHomCols)
#     kinship_dev_vec <- abs((2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]] + SUPER_15_peak1_GL[pc1_filter,MajorHomCols[i]+1]) - (2*SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]]+ SUPER_15_peak1_GL[pc1_filter,MajorHomCols[1]+1]))
#     kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]] + mean(kinship_avg_list[[Beagle_ind_list[i,"pc1_clust"]]], na.rm = T)
#     #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
#     if(i != 1) lines(x = SUPER_15_peak1_GL[pc1_filter,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3)+1, lwd = 0.6, col = rgb(t(col_vec_hl[,Beagle_ind_list[i,"pc1_clust"]]), maxColorValue = 255))
#     #Individual 1 is the comparison reference.
#   }
# }
# dev.off()

# #Frequency plots
# pdf("~/Projects/Eel/doc/SUPER_15_peak1_freq.pdf", width = 10, height = 10)
# plot(x= all_freq_avg_list[[1]], y = all_freq_avg_list[[2]], col = "grey70", cex = 0.5, xlab = "Ref. allele freq. in T1 homozygotes", ylab = "Ref. allele freq. in T2 homozygotes", main = "SUPER_15: 12.1 Mb")
# shared_snps <- abs(all_freq_avg_list[[1]] - all_freq_avg_list[[2]]) < 0.2 & all_freq_avg_list[[1]] < 0.75 & all_freq_avg_list[[2]] < 0.75
# points(x = all_freq_avg_list[[1]][shared_snps], y = all_freq_avg_list[[2]][shared_snps], col = "darkorchid", pch = 20, cex = 0.5)
# points(x = pc1_high_freq_avg_list[[1]], y = pc1_high_freq_avg_list[[2]], col = "red", pch = 20, cex = 0.5)
# legend(x = "topleft", legend = c("Shared", "Diagnostic", "Background"), col = c("darkorchid", "red", "grey70"), pch = c(20, 20, 1))
# 
# plot(y= abs(all_freq_avg_list[[1]] - all_freq_avg_list[[2]]), x = SUPER_15_peak1_GL[,"pos"], col = "grey70", cex = 0.5, xlab = "Position", main = "SUPER_15", ylab = "DAF (T1 hom vs T2 hom)")
# points(x = SUPER_15_peak1_GL[,"pos"][shared_snps], y = abs(all_freq_avg_list[[1]] - all_freq_avg_list[[2]])[shared_snps], col = "darkorchid", pch = 20, cex = 0.5)
# points(x = SUPER_15_peak1_GL[pc1_filter,"pos"], y = abs(pc1_high_freq_avg_list[[2]] - pc1_high_freq_avg_list[[1]]), col = "red", pch = 20, cex = 0.5)
# legend(x = "topleft", legend = c("Shared", "Diagnostic", "Background"), col = c("darkorchid", "red", "grey70"), pch = c(20, 20, 1))
# dev.off()



#Support functions
kin_and_freq <- function(reg_GL, clust_vec,  ref_ind = 1){ #Adapted to deal with NAs in the GL data #pdf_file, old argument?
  MajorHomCols <- grep("Ind[0-9]+$", names(reg_GL))
  
  #Marker-wise stat collection, allowing for arbitrary numbers of individual classes
  n_groups <- max(clust_vec)
  gl_avg_df <- data.frame(type1_kin = numeric(dim(reg_GL)[1]))
  for(c in 1:n_groups){
    gl_avg_df[,paste0("type", c, "_kin")] <- numeric(dim(reg_GL)[1])
  }
  for(c in 1:n_groups){
    gl_avg_df[,paste0("type", c, "_freq")] <- numeric(dim(reg_GL)[1])
  }
  for(c in 1:n_groups){
    gl_avg_df[,paste0("type", c, "_nobs")] <- numeric(dim(reg_GL)[1])
  }
                          
  #gl_avg_df <- data.frame(t1_hom_kin = numeric(dim(reg_GL)[1]), t2_hom_kin = numeric(dim(reg_GL)[1]), het_kin  = numeric(dim(reg_GL)[1]), t1_hom_freq = numeric(dim(reg_GL)[1]), t2_hom_freq = numeric(dim(reg_GL)[1]), het_freq = numeric(dim(reg_GL)[1]), t1_hom_nobs = numeric(dim(reg_GL)[1]), t2_hom_nobs = numeric(dim(reg_GL)[1]), het_nobs = numeric(dim(reg_GL)[1]))
  rownames(gl_avg_df) <- reg_GL$marker
  #all_freq_avg_df <- list(t1_hom = numeric(dim(reg_GL)[1]), t2_hom = numeric(dim(reg_GL)[1]), het = numeric(dim(reg_GL)[1]))
  
  for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
    kinship_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
    freq_vec <- abs(2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1])
    gl_avg_df[which(!is.na(freq_vec)),clust_vec[i]+n_groups*2] <- gl_avg_df[which(!is.na(freq_vec)),clust_vec[i]+n_groups*2] + 1 #collect no of observations
    kinship_vec[is.na(kinship_vec)] <- 0  #The set NAs to zero, to enable summation
    freq_vec[is.na(freq_vec)] <- 0  #The set NAs to zero, to enable summation
    gl_avg_df[,clust_vec[i]] <- gl_avg_df[,clust_vec[i]] +  kinship_vec #kinship
    gl_avg_df[,clust_vec[i]+n_groups] <- gl_avg_df[,clust_vec[i]+n_groups] + freq_vec #greq
    #all_freq_avg_list[[clust_vec[i]]] <- all_freq_avg_list[[clust_vec[i]]] + (2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1])
  }
  
  for(i in 1:n_groups){
    gl_avg_df[gl_avg_df[,i+n_groups*2] == 0, i+n_groups*2] <- NA #Set non-observed classes to NA
    gl_avg_df[,i] <- gl_avg_df[,i]/gl_avg_df[,i+2*n_groups]
    gl_avg_df[,i+n_groups] <- gl_avg_df[,i+n_groups]/(2*gl_avg_df[,i+n_groups*2])
   #Kinship is onlt meaningful if reference individual is observed
    ref_filter <- is.na(reg_GL[,MajorHomCols[ref_ind]])
    gl_avg_df[ref_filter,i] <- NA
    #kinship_avg_list[[i]] <- kinship_avg_list[[i]]/sum(clust_vec == i, na.rm = T)
    #all_freq_avg_list[[i]] <- all_freq_avg_list[[i]]/(2*sum(clust_vec == i, na.rm = T))
  }
  
  
  #Needs update to new format
  # col_vec <- col2rgb(c("blue", "red",  "darkorchid"),alpha = T)
  # col_vec["alpha",] <- 80
  # 
  # pdf(pdf_file, width = 10)
  # plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "", xlab = "Position")
  # for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
  #   kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
  #   kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
  #   #Adjusting for the group avearge at each marker, and also offsetting according to the mean
  #   if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.3, col = rgb(t(col_vec[,clust_vec[i]]), maxColorValue = 255)) 
  # }
  # 
  # col_vec_hl <- col_vec
  # col_vec_hl["alpha",] <- 150
  # 
  # for(j in 1:3){
  #   plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "", xlab = "Position")
  #   for(i in which(clust_vec != j)){ #length(MajorHomCols)
  #     kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
  #     kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
  #     #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
  #     if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.2, col = "grey70")
  #     #Individual 1 is the comparison reference.
  #   }
  #   for(i in which(clust_vec == j)){ #length(MajorHomCols)
  #     kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
  #     kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
  #     #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
  #     if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3)+1, lwd = 0.6, col = rgb(t(col_vec_hl[,clust_vec[i]]), maxColorValue = 255))
  #     #Individual 1 is the comparison reference.
  #   }
  # }
  # dev.off()
  return(invisible(gl_avg_df))
}

sliding_deviation <- function(reg_GL, clust_vec, pdf_file, ref_ind = 1){
  MajorHomCols <- grep("Ind[0-9]+$", names(reg_GL))
  
  kinship_avg_list <- list(t1_hom = numeric(dim(reg_GL)[1]), t2_hom = numeric(dim(reg_GL)[1]), het = numeric(dim(reg_GL)[1]))
  all_freq_avg_list <- list(t1_hom = numeric(dim(reg_GL)[1]), t2_hom = numeric(dim(reg_GL)[1]), het = numeric(dim(reg_GL)[1]))
  obs_ind_df <- data
  for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
    kinship_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
    kinship_avg_list[[clust_vec[i]]] <- kinship_avg_list[[clust_vec[i]]] + kinship_vec
    all_freq_avg_list[[clust_vec[i]]] <- all_freq_avg_list[[clust_vec[i]]] + (2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1])
  }
  for(i in 1:3){
    kinship_avg_list[[i]] <- kinship_avg_list[[i]]/sum(clust_vec == i, na.rm = T)
    all_freq_avg_list[[i]] <- all_freq_avg_list[[i]]/(2*sum(clust_vec == i, na.rm = T))
  }
  
  
  
  col_vec <- col2rgb(c("blue", "red",  "darkorchid"),alpha = T)
  col_vec["alpha",] <- 80
  
  pdf(pdf_file, width = 10)
  plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "", xlab = "Position")
  for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
    kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
    kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
    #Adjusting for the group avearge at each marker, and also offsetting according to the mean
    if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.3, col = rgb(t(col_vec[,clust_vec[i]]), maxColorValue = 255)) 
  }
  
  col_vec_hl <- col_vec
  col_vec_hl["alpha",] <- 150
  
  for(j in 1:3){
    plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = "Deviaiton from group expectation (100 SNP rolling avg)", main = "", xlab = "Position")
    for(i in which(clust_vec != j)){ #length(MajorHomCols)
      kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
      kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
      #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
      if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3) +1, lwd = 0.2, col = "grey70")
      #Individual 1 is the comparison reference.
    }
    for(i in which(clust_vec == j)){ #length(MajorHomCols)
      kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[ref_ind]]+ reg_GL[,MajorHomCols[ref_ind]+1]))
      kinship_dev_vec <- kinship_dev_vec - kinship_avg_list[[clust_vec[i]]] + mean(kinship_avg_list[[clust_vec[i]]], na.rm = T)
      #Adjusting for the group avearge at each marker, and also offsetting the groups according to their mean
      if(i != ref_ind) lines(x = reg_GL[,"pos"], y = filter(kinship_dev_vec, rep(1/100, 100)), lty = (i %% 3)+1, lwd = 0.6, col = rgb(t(col_vec_hl[,clust_vec[i]]), maxColorValue = 255))
      #Individual 1 is the comparison reference.
    }
  }
  dev.off()
  return(invisible(list(kin = kinship_avg_list, freq = all_freq_avg_list)))
}

region_kinship <- function(reg_GL, ind_list){ #
  MajorHomCols <- grep("Ind[0-9]+$", names(reg_GL))
  kinship_mat <- matrix(nrow = length(ind_list[,"ID"]), ncol = length(ind_list[,"ID"]))
  rownames(kinship_mat) <- ind_list[,"ID"]
  colnames(kinship_mat) <- ind_list[,"ID"]
  #kinship_mat_filter <- NULL
  #if(!is.null(marker_filter[1])){
  #  kinship_mat_high_pc1 <- matrix(nrow = length(ind_list[,"ID"]), ncol = length(ind_list[,"ID"]))
  #  rownames(kinship_mat_filter) <- ind_list[,"ID"]
  #  colnames(kinship_mat_filter) <- ind_list[,"ID"]
  #}
  
  for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
    for(j in i:length(MajorHomCols)){ #length(MajorHomCols)
      kinship_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[j]] + reg_GL[,MajorHomCols[j]+1]))
      kinship_mat[i,j] <- sum(kinship_vec, na.rm = T)/sum(!is.na(kinship_vec))
      kinship_mat[j,i] <- kinship_mat[i,j]
      #if(!is.null(marker_filter[1])){
      #  kinship_mat_filter[i,j] <- sum(abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - (2*reg_GL[,MajorHomCols[j]] + reg_GL[,MajorHomCols[j]+1]))[marker_filter])/sum(marker_filter)
      #  kinship_mat_filter[j,i] <- kinship_mat_filter[i,j]
      #}
    }
  }
  return(list(full=kinship_mat)) # kinship_mat_filter
}


plot_region_hm <- function(target_GR, g1_ids = NULL, g2_ids = NULL, ind_co = 350, daf_co = 0.15, plot_dir = "~/Projects/Eel/doc/", target_freq = NULL){
  t_super <- as.character(target_GR@seqnames)
  #if(is.null(target_freq)){
  #  target_freq<- get(paste0(t_super,"_freq"))
  #}
  pos_filter <- target_freq$nInd > ind_co & target_freq$position > target_GR@ranges@start & target_freq$position < target_GR@ranges@start + target_GR@ranges@width
  t_reg_freq <- target_freq[pos_filter,-(1:8)]
  rownames(t_reg_freq) <- target_freq[pos_filter,"position"]
  high_daf <- !logical(length = dim(t_reg_freq)[1])
  pdf_file <- gsub("[:-]", "_", paste0(plot_dir, paste(target_GR), "_hm.pdf"))
  pdf(pdf_file)
  if(!is.null(g1_ids) & !is.null(g2_ids)){
    t_daf <- abs(rowMeans(as.data.frame(t_reg_freq[,g1_ids]), na.rm = T) - rowMeans(t_reg_freq[,g2_ids], na.rm = T))
    high_daf <- t_daf > daf_co
    plot(y = t_daf, x = rownames(t_reg_freq), xlab = "Position", main = paste0(target_GR@seqnames[1], ": ", paste(g1_ids,collapse = "+"), " vs ", paste(g2_ids,collapse = "+")), ylab = "DAF", pch = 20, cex = 0.8)
    abline(h = daf_co, lty = "dashed", col = "red")
    heatmap(t(as.matrix(t_reg_freq[high_daf,])), scale = "none", Colv = NA, main = paste(paste(target_GR), "High DAF SNPs"))
  }
  heatmap(t(as.matrix(t_reg_freq)), scale = "none", Colv = NA, main = paste(paste(target_GR), "All SNPs"))
  dev.off()
  
}

plot_contrast <- function(g1_ids, g2_ids, contrast_folder = "g1_v_g2", ind_co = 350, chr_list = 1:19, var_suffix = "_0.05_freq"){
  #total_out_df <- NULL
  chr_vec <- character()
  pos_vec <- integer()
  cpos_vec <- integer()
  r_DAF_vec <- numeric()
  DAF_vec <- numeric()
  col_vec <- character()
  pos_offset <- 0
  col_counter <- 1
  col_set <- c("grey30", "grey70")
  for(t_super in paste0("SUPER_", chr_list)){
    #if(!exists(paste0(t_super,"_freq"))){
    if(!exists(paste0(t_super, var_suffix))){
      #load(paste0("~/Projects/Eel/data/all_maf/",t_super,"_freq.RData"))
      load(paste0("~/Projects/Eel/data/all_maf/maf_co_0.05/",t_super,var_suffix,".RData"))
    }
    #target_freq <- get(paste0(t_super,"_freq"))
    target_freq <- get(paste0(t_super,var_suffix))
    print(paste("Target variable:", paste0(t_super,var_suffix)))
    flush.console()
    pos_filter <- target_freq$nInd > ind_co
    t_daf <- abs(rowMeans(as.data.frame(target_freq[,g1_ids]), na.rm = T) - rowMeans(target_freq[,g2_ids], na.rm = T))
    #if(is.null(x_lim)){
    #png(paste0("~/Projects/Eel/doc/",contrast_folder,"/", t_super,".png"), width = 1000)
    png(paste0(contrast_folder,"/", t_super,".png"), width = 1000)
    plot(x = target_freq$position[pos_filter], y =  t_daf[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = t_super, col = "grey70", ylim=c(0,1))
    #}
    
    #if(!is.null(x_lim)){
    #  target_reg <- round(x_lim/1e6, digits = 1)
    #  png(paste0("~/Projects/Eel/doc/SWE_v_EU/", target_chr,"_", target_reg[1], "_to_", target_reg[2], "Mb_SWEvEU.png"), width = 1000)
    #  plot(x = target_chr_maf$position[pos_filter], y = target_chr_maf$DAF[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = target_chr, col = "grey70", ylim=c(0,1), xlim = x_lim)
    # 
    #}
    #plot(x = target_chr_maf$position[pos_filter], y = target_chr_maf$DAF[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = target_chr, col = "grey70", ylim=c(0,1))
    rolling_DAF <- stats::filter(x= t_daf[pos_filter], rep(1/100, 100))
    points(x = target_freq$position[pos_filter], y = rolling_DAF, pch = 16, cex = 0.5, ylab = "DAF (100 SNP average)", xlab = "Position", col = "firebrick")
    dev.off()
    
    #if(t_super == )){
      #total_out_df <- data.frame(chr = t_super, DAF = t_daf[pos_filter], r_DAF = rolling_DAF, pos = target_freq$position[pos_filter], col = col_vec[(col_counter %% 2) +1], stringsAsFactors = F)
    #  col_counter <- col_counter + 1
    #} else{
      #tmp_df <- data.frame(chr = t_super, DAF = t_daf[pos_filter], r_DAF = rolling_DAF, pos = target_freq$position[pos_filter], col = col_vec[(col_counter %% 2) +1], stringsAsFactors = F)
      #total_out_df <- rbind(total_out_df, tmp_df)
    chr_vec <- c(chr_vec, rep(t_super, sum(pos_filter)))
    col_vec <- c(col_vec, rep(col_set[(col_counter %% 2) +1], sum(pos_filter)))
    DAF_vec <- c(DAF_vec, t_daf[pos_filter])
    r_DAF_vec <- c(r_DAF_vec, rolling_DAF)
    pos_vec <- c(pos_vec, target_freq$position[pos_filter])
    cpos_vec <- c(cpos_vec, (target_freq$position[pos_filter] + pos_offset))
    
    pos_offset <- pos_offset + max(target_freq$position[pos_filter])           
    col_counter <- col_counter + 1
    rm(list = paste0(t_super,var_suffix))
    #}
  }
 
  png(paste0(contrast_folder,"/GW_DAF.png"), width = 10000, height = 1000)
  plot(x = cpos_vec, y =  DAF_vec, pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = "", col = col_vec, ylim=c(0,1))
  points(x = cpos_vec, y = r_DAF_vec, pch = 16, cex = 0.5, ylab = "DAF (100 SNP average)", xlab = "Position", col = "firebrick")
  dev.off()
  
  #return(invisible(total_out_df))
  return(invisible(list(chr = chr_vec, pos = pos_vec, cpos = cpos_vec, DAF = DAF_vec, rDAF = r_DAF_vec, col = col_vec)))
}


read_maf_files <- function(target_chr, group_ids, total_dir = "", per_group_dir = ""){
  
  #ALL_file <- paste0("~/Projects/Eel/data/all_maf/All_MnM_", target_chr, ".mafs")
  ALL_file <- grep(paste0(target_chr,"[._]maf"), dir(total_dir, full.names = T), value = T)
  target_chr_maf <- read.table(ALL_file, header = T, stringsAsFactors = F, sep = "\t")
  target_chr_maf[,"SNP_ID"] <- paste(target_chr_maf$chromo,target_chr_maf$position, sep = "_")
  for(g_name in group_ids){
    #group_file <- paste0("~/Projects/Eel/data/per_group_maf/", g_name, "_MnM_", target_chr, ".mafs.gz")
    group_file <- grep(paste0(g_name,".+", target_chr, "_.+mafs[.]gz"),dir(per_group_dir, full.names = T), value = T)
    maf_group <- read.table(group_file, header = T, stringsAsFactors = F, sep = "\t")
    maf_group[,"SNP_ID"] <- paste(maf_group$chromo, maf_group$position, sep = "_")
    target_chr_maf[,g_name] <- maf_group[match(target_chr_maf$SNP_ID, maf_group$SNP_ID), "knownEM"]
  }
  return(invisible(target_chr_maf))
}

read_and_plot_contrast <- function(target_chr, ind_co = 400, x_lim = NULL){
  ALL_file <- paste0("~/Projects/Eel/data/SWE_v_EU_maf/All_MnM_", target_chr, ".mafs")
  SWE_file <- paste0("~/Projects/Eel/data/SWE_v_EU_maf/SWE_MnM_", target_chr, ".mafs.gz")
  EU_file <- paste0("~/Projects/Eel/data/SWE_v_EU_maf/EU_MnM_", target_chr, ".mafs.gz")
  
  target_chr_maf <- read.table(ALL_file, header = T, stringsAsFactors = F, sep = "\t")
  target_chr_maf[,"SNP_ID"] <- paste(target_chr_maf$chromo,target_chr_maf$position, sep = "_")
  
  maf_SWE <- read.table(SWE_file, header = T, stringsAsFactors = F, sep = "\t")
  maf_SWE[,"SNP_ID"] <- paste(maf_SWE$chromo,maf_SWE$position, sep = "_")
  
  maf_EU <- read.table(EU_file, header = T, stringsAsFactors = F, sep = "\t")
  maf_EU[,"SNP_ID"] <- paste(maf_EU$chromo,maf_EU$position, sep = "_")
  
  target_chr_maf[,c("SWE_maf", "SWE_ID")] <- maf_SWE[match(target_chr_maf$SNP_ID, maf_SWE$SNP_ID), c("knownEM", "SNP_ID")]
  target_chr_maf[,c("EU_maf", "EU_ID")] <- maf_EU[match(target_chr_maf$SNP_ID, maf_EU$SNP_ID), c("knownEM", "SNP_ID")] 
  
  target_chr_complete <- !is.na(target_chr_maf$SWE_maf) & !is.na(target_chr_maf$EU_maf)
  
  target_chr_maf <- target_chr_maf[target_chr_complete,]
  #sum(SUPER_19_maf$SNP_ID == SUPER_19_maf$EU_ID & SUPER_19_maf$SNP_ID == SUPER_19_maf$SWE_ID) #Checks out
  target_chr_maf[,"DAF"] <- abs(target_chr_maf$SWE_maf - target_chr_maf$EU_maf)
  
  #Visualisation
  #ind_co <- 400
  pos_filter <- target_chr_maf$nInd > ind_co
  if(is.null(x_lim)){
    png(paste0("~/Projects/Eel/doc/SWE_v_EU/", target_chr, "_SWEvEU.png"), width = 1000)
    plot(x = target_chr_maf$position[pos_filter], y = target_chr_maf$DAF[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = target_chr, col = "grey70", ylim=c(0,1))
  }
  if(!is.null(x_lim)){
    target_reg <- round(x_lim/1e6, digits = 1)
    png(paste0("~/Projects/Eel/doc/SWE_v_EU/", target_chr,"_", target_reg[1], "_to_", target_reg[2], "Mb_SWEvEU.png"), width = 1000)
    plot(x = target_chr_maf$position[pos_filter], y = target_chr_maf$DAF[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = target_chr, col = "grey70", ylim=c(0,1), xlim = x_lim)
    
  }
  #plot(x = target_chr_maf$position[pos_filter], y = target_chr_maf$DAF[pos_filter], pch = 16, cex = 0.5, ylab = "DAF", xlab = "Position", main = target_chr, col = "grey70", ylim=c(0,1))
  rolling_DAF <- filter(x=target_chr_maf$DAF[pos_filter], rep(1/100, 100))
  points(x = target_chr_maf$position[pos_filter], y = rolling_DAF, pch = 16, cex = 0.5, ylab = "DAF (100 SNP average)", xlab = "Position", col = "firebrick")
  dev.off()
  return(invisible(target_chr_maf ))
}
