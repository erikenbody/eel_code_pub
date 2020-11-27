#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Regions of interest

#SUPER_15 ~12 Mb
SUPER_15_site_sel <- read.csv(file = "~/Projects/Eel/data/geno_likelihood/selection_SUPER_15.csv", stringsAsFactors = F)
SUPER_15_site_sel[,"position"] <- as.numeric(sub("SUPER_15_", "", SUPER_15_site_sel[,"pos"]))
#SUPER 15  peak 1
#SUPER_15_11870507 - SUPER_15_12246140
#zgrep -m1 -n -o SUPER_15_11870507 SUPER_19_SUPER_15.beagle.gz
325503:SUPER_15_11870507
332246:SUPER_15_12246140
#zcat SUPER_19_SUPER_15.beagle.gz | awk 'NR == 1 || (NR >=325503 && NR <= 332246) {print}' > SUPER_15_peak1.beagle 

SUPER_15_peak1_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak1.beagle.gz", header = T, sep = "\t", stringsAsFactors = F)
SUPER_15_peak1_GL[,"pos"] <- as.numeric(sub("SUPER_15_", "", SUPER_15_peak1_GL[,"marker"]))
SUPER_15_peak1_GL[,"pc1.pval"] <- SUPER_15_site_sel[match(SUPER_15_peak1_GL[,"pos"], SUPER_15_site_sel[,"position"]),"pc1.pval"]

SUPER_15_GL_clean <- SUPER_15_peak1_GL
MajorHomCols <- grep("Ind[0-9]+$", names(SUPER_15_peak1_GL))
pc1_filter <- SUPER_15_peak1_GL$pc1.pval < 1e-10
pc1_filter[is.na(pc1_filter)] <- F
for(i in 1:length(MajorHomCols)){
  ind_SDs <- rowSds(as.matrix(SUPER_15_GL_clean[, MajorHomCols[i] + c(0,1,2)]), na.rm = T)
  SUPER_15_GL_clean[ind_SDs == 0,  MajorHomCols[i] + c(0,1,2)] <- NA
}

SUPER_15_kinship <- region_kinship(SUPER_15_GL_clean[pc1_filter,], Beagle_ind_list)
SUPER_15_kinship[[1]][is.na(SUPER_15_kinship[[1]])] <- mean(SUPER_15_kinship[[1]], na.rm = T) #Some pairs APM45/AMAR10 have no value, setting them to the overall avearge
SUPER_15_hclust <- hclust(as.dist(SUPER_15_kinship[[1]]))
plot(SUPER_15_hclust)
SUPER_15_c_assign <- rect.hclust(SUPER_15_hclust, k = 4)
for(i in 1:4) Beagle_ind_list[SUPER_15_c_assign[[i]],"super_15_clust_v2"] <- i
table(Beagle_ind_list[,c("pc1_clust", "super_15_clust_v2")])

save(Beagle_ind_list, file = "~/Projects/Eel/data/individual_group_assignment.RData")

SUPER_15_freq <- kin_and_freq(SUPER_15_GL_clean, Beagle_ind_list$super_15_clust_v2)
save(SUPER_15_freq, SUPER_15_kinship, SUPER_15_GL_clean, file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak1_clean.RData")

pdf(file = "~/Projects/Eel/doc/SUPER_15_freq_update.pdf")
heatmap(SUPER_15_kinship[[1]], scale = "none", Rowv = as.dendrogram(SUPER_15_hclust), as.dendrogram(SUPER_15_hclust))
plot(x = SUPER_15_freq$type4_freq[!pc1_filter], y = SUPER_15_freq$type1_freq[!pc1_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3, xlab = "Frequency in Type 1 homozygotes",  ylab = "Frequency in Type 2 homozygotes") 
points(x = SUPER_15_freq$type4_freq[pc1_filter], y = SUPER_15_freq$type1_freq[pc1_filter], xlim = c(0,1), ylim = c(0,1), pch = 16, col = "darkorchid", cex = 0.3)

plot(x = SUPER_15_GL_clean$pos[!pc1_filter], y = abs(SUPER_15_freq$type4_freq[!pc1_filter]-SUPER_15_freq$type1_freq[!pc1_filter]), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_15_GL_clean$pos[pc1_filter], y = abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]), pch = 16, col = "darkorchid", cex = 0.3)
dev.off()

#Output high DAF markers to intersect with snpEFF file
write(SUPER_15_GL_clean$pos[abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]) > 0.8], file = "~/Projects/Eel/data/annotation/SUPER_15_high_DAF.txt", ncolumns = 1)

#Diverisity calculations - probabilities of different alleles  
p_1ne4 <- 1-(SUPER_15_freq$type1_freq * SUPER_15_freq$type4_freq + (1-SUPER_15_freq$type1_freq) * (1-SUPER_15_freq$type4_freq))
p_1ne1 <- 1-(SUPER_15_freq$type1_freq * SUPER_15_freq$type1_freq + (1-SUPER_15_freq$type1_freq) * (1-SUPER_15_freq$type1_freq))
p_4ne4 <- 1-(SUPER_15_freq$type4_freq * SUPER_15_freq$type4_freq + (1-SUPER_15_freq$type4_freq) * (1-SUPER_15_freq$type4_freq))
#Call rate filter
high_call <- rowSums(is.na(SUPER_15_GL_clean[,MajorHomCols])) < 150

rolling_4ne4 <- filter(p_4ne4, rep(1, 100))
rolling_1ne1 <- filter(p_1ne1, rep(1, 100))
rolling_1ne4 <- filter(p_1ne4, rep(1, 100))
rolling_pos <- c(rep(NA, 49), diff(SUPER_15_GL_clean$pos, lag = 99), rep(NA, 50))

pdf(file = "~/Projects/Eel/doc/SUPER_15_divergence.pdf", width = 10)
plot(x = SUPER_15_GL_clean$pos[high_call], y = p_1ne4[high_call], xlab = "Position", ylab = "Probability of finding different allele", pch = 16, col = "grey20", ylim = c(0,1), cex = 0.3, main = "SUPER_15")
points(x = SUPER_15_GL_clean$pos[high_call], y = p_4ne4[high_call], pch = 16, col = "blue", cex = 0.3)
legend(x="topleft", legend = c("Type 1 vs Type 2","Within Type 1"), col = c("black", "blue"), pch = 16)

plot(x = SUPER_15_GL_clean$pos[high_call], y = p_1ne4[high_call], xlab = "Position", ylab = "Probability of finding different allele", pch = 16, col = "grey20", ylim = c(0,1), cex = 0.3, main = "SUPER_15")
points(x = SUPER_15_GL_clean$pos[high_call], y = p_1ne1[high_call], pch = 16, col = "red", cex = 0.3)
legend(x="topleft", legend = c("Type 1 vs Type 2","Within Type 2"), col = c("black", "red"), pch = 16)

plot(x = SUPER_15_GL_clean$pos, y = (rolling_1ne4/rolling_pos), xlab = "Position", ylab = "Diversity (100 SNP rolling avg)", pch = 16, col = "grey20", lwd = 0.7, type = "l", main = "SUPER_15")
lines(x = SUPER_15_GL_clean$pos, y = (rolling_1ne1/rolling_pos), pch = 16, col = "red", lwd = 0.5)
lines(x = SUPER_15_GL_clean$pos, y = (rolling_4ne4/rolling_pos), pch = 16, col = "blue", lwd = 0.5)
legend(x="topright", legend = c("Type 1 vs Type 2","Within Type 1", "Within Type 2"), col = c("black", "blue", "red"), lty = 1)
dev.off()


sum(p_4ne4)/diff(range(SUPER_15_GL_clean$pos))
sum(p_1ne4, na.rm = T)/diff(range(SUPER_15_GL_clean$pos))
sum(p_1ne1, na.rm = T)/diff(range(SUPER_15_GL_clean$pos))

#Type frequencies per group
SUPER_15_p1_by_group <- table(Beagle_ind_list$super_15_clust_v2, Beagle_ind_list$GROUP)
(SUPER_15_p1_by_group[1,] + SUPER_15_p1_by_group[2,]*0.5)/colSums(SUPER_15_p1_by_group[-3,])
#AAF       AES        AI       ALI      AMAR       APM      ASRI       AST       ATU       SMS 
#0.3750000 0.2678571 0.2708333 0.2291667 0.3125000 0.2500000 0.2666667 0.2833333 0.2761194 0.3723404
#sum(SUPER_15_p1_by_group[1,] + SUPER_15_p1_by_group[2,]*0.5)/445 #overall mean
#[1] 0.2898876
exp_vec <- 0.2898876*colSums(SUPER_15_p1_by_group[,])*2
obs_vec <- (2*SUPER_15_p1_by_group[1,] + SUPER_15_p1_by_group[2,])
test_mat <- matrix(data = c(obs_vec, exp_vec), ncol = 2)
chisq.test(test_mat)

#Annotation subset & sliding plot
SUPER_15_GR <- GRanges(seqnames = (fAngAng_gff[fAngAng_gff$type == "region"][15])@seqnames, ranges = IRanges(start = 11.8e6, end = 12.2e6))
SUPER_15_hits <- findOverlaps(SUPER_15_GR, fAngAng_gff)
SUPER_15_gff <- fAngAng_gff[SUPER_15_hits@to]
SUPER_15_gff[SUPER_15_gff$type == "gene"]

#sliding_deviation_v2(reg_GL = SUPER_15_GL_clean, ref_freq_vec = SUPER_15_freq$type1_freq, clust_vec = Beagle_ind_list$super_15_clust_v2, pdf_file = "~/Projects/Eel/doc/SUPER_15_p1_slide_update.pdf")
sliding_deviation_v2(reg_GL = SUPER_15_GL_clean, ref_freq_vec = SUPER_15_freq$type1_freq, clust_vec = Beagle_ind_list$super_15_clust_v2, png_file = "~/Projects/Eel/doc/SUPER_15_p1_slide_update.png", target_gff = SUPER_15_gff)
sliding_deviation_v2(reg_GL = SUPER_15_GL_clean[pc1_filter,], ref_freq_vec = SUPER_15_freq$type1_freq[pc1_filter], clust_vec = Beagle_ind_list$super_15_clust_v2, pdf_file = "~/Projects/Eel/doc/SUPER_15_p1_slide_diag.pdf", smooth_int = 25)

#Missense mutations
SUPER_15_snpEff <- read.table("~/Projects/Eel/data/annotation/SUPER_15_high_DAF_snpEff.out", stringsAsFactors = F, sep = "\t")
SUPER_15_snpEff_high_DAF <- SUPER_15_snpEff[SUPER_15_snpEff$V2 %in% SUPER_15_GL_clean$pos[pc1_filter],]
#Alternative set, based on DAF between types
S15_p1_daf_filter <- abs(SUPER_15_freq$type1_freq - SUPER_15_freq$type4_freq) > 0.6 & SUPER_15_freq$type1_nobs > 20 &  SUPER_15_freq$type4_nobs > 125
SUPER_15_snpEff_high_DAF_alt <- SUPER_15_snpEff[SUPER_15_snpEff$V2 %in% SUPER_15_GL_clean$pos[S15_p1_daf_filter],]
S15_miss_rows <- grep("missense", SUPER_15_snpEff_high_DAF_alt$V8)
S15_snpEff_high_DAF_alt_miss <- SUPER_15_snpEff_high_DAF_alt[S15_miss_rows,]

SUPER_15_GR <- GRanges(seqnames = (fAngAng_gff[fAngAng_gff$type == "region"][15])@seqnames, ranges = IRanges(start = 11.8e6, end = 12.2e6))
SUPER_15_hits <- findOverlaps(SUPER_15_GR, fAngAng_gff)
SUPER_15_gff <- fAngAng_gff[SUPER_15_hits@to]
SUPER_15_gff[SUPER_15_gff$type == "gene"]

SUPER_15_missense_AA <- as.integer(sub(".+p[.][A-Za-z]{3}([0-9]+)[A-Za-z]{3}.+", "\\1", SUPER_15_snpEff_high_DAF[grep("missense", SUPER_15_snpEff_high_DAF$V8),8]))
names(SUPER_15_missense_AA) <- sub(".+rna[-](XM_[0-9]+[.][0-9]).+", "\\1", SUPER_15_snpEff_high_DAF[grep("missense", SUPER_15_snpEff_high_DAF$V8),8])
table(names(SUPER_15_missense_AA))

pdf(file = "~/Projects/Eel/doc/SUPER_15_clean_annot.pdf", width = 10)
plot(x = SUPER_15_GL_clean$pos[!pc1_filter], y = abs(SUPER_15_freq$type4_freq[!pc1_filter]-SUPER_15_freq$type1_freq[!pc1_filter]), ylim = c(0,1), pch = 16, col = "grey70", cex = 0.3, ylab = "DAF", xlab = "Position")
points(x = SUPER_15_GL_clean$pos[pc1_filter], y = abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]), pch = 16, col = "cornflowerblue", cex = 0.4)
#abline(v = SUPER_15_snpEff_high_DAF[grep("missense", SUPER_15_snpEff_high_DAF$V8),2])
segments(x0 = S15_snpEff_high_DAF_alt_miss[,2], y0 = abs(SUPER_15_freq$type4_freq[S15_p1_daf_filter]-SUPER_15_freq$type1_freq[S15_p1_daf_filter])[S15_miss_rows], y1 = 0.9, col = "grey70", lwd = 1, lty = "dashed")
points(x= S15_snpEff_high_DAF_alt_miss[,2], y = abs(SUPER_15_freq$type4_freq[S15_p1_daf_filter]-SUPER_15_freq$type1_freq[S15_p1_daf_filter])[S15_miss_rows], pch = 20, col = "red")

segments(x0 = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start, x1 = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start+SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@width, y0 = 0.9, col = "firebrick", lwd = 3)
#abline(v = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start,  col = "firebrick1", lty = "dashed", lwd = 0.5)
#abline(v = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start + SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@width, col = "firebrick1", lty = "dashed", lwd = 0.5)

x0_vec <- SUPER_15_gff[SUPER_15_gff$type == "exon"]@ranges@start
x1_vec <- x0_vec + SUPER_15_gff[SUPER_15_gff$type == "exon"]@ranges@width
rect(xleft = x0_vec, xright = x1_vec, ybottom = 0.88, ytop = 0.92, col = "firebrick", border = NA)

gene_mid_vec <- SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start+ (0.5*SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@width)
text(x = gene_mid_vec[2:5], y = c(0.95, 0.85), labels = SUPER_15_gff[SUPER_15_gff$type == "gene"]$Name[2:5], cex = 1)

dev.off()

#All high-DAF mutations from higly affected transcript "XM_035392574.1"
XM_035392574_mut <- SUPER_15_snpEff_high_DAF_alt[grep("XM_035392574[.]1", SUPER_15_snpEff_high_DAF_alt$V8), ]
sum(grepl("syno", XM_035392574_mut$V8))
#[1] 13
sum(grepl("missense", XM_035392574_mut$V8))
#[1] 17
#"XM_035394417.1" 
XM_035394417_mut <- SUPER_15_snpEff_high_DAF_alt[grep("XM_035394417[.]1", SUPER_15_snpEff_high_DAF_alt$V8), ]
sum(grepl("syno", XM_035394417_mut$V8))
#[1] 8
sum(grepl("missense", XM_035394417_mut$V8))
#[1] 1

#"XM_035393870.1"
XM_035393870_mut <- SUPER_15_snpEff_high_DAF_alt[grep("XM_035393870[.]1", SUPER_15_snpEff_high_DAF_alt$V8), ]
sum(grepl("syno", XM_035393870_mut$V8))
#[1] 15
sum(grepl("missense", XM_035393870_mut$V8))
#[1] 3



tmp_S15_GL <- SUPER_15_GL_clean[S15_p1_daf_filter,1:10]
tmp_S15_freq <- SUPER_15_0.05_freq[SUPER_15_0.05_freq$position %in% SUPER_15_GL_clean$pos[S15_p1_daf_filter],]
tmp_S15_GL[,"a1_letter"] <- "N"
tmp_S15_GL[tmp_S15_GL$allele1 == 0,"a1_letter"] <- "A"
tmp_S15_GL[tmp_S15_GL$allele1 == 1,"a1_letter"] <- "C"
tmp_S15_GL[tmp_S15_GL$allele1 == 2,"a1_letter"] <- "G"
tmp_S15_GL[tmp_S15_GL$allele1 == 3,"a1_letter"] <- "T"
tmp_S15_GL[,"a2_letter"] <- "N"
tmp_S15_GL[tmp_S15_GL$allele2 == 0,"a2_letter"] <- "A"
tmp_S15_GL[tmp_S15_GL$allele2 == 1,"a2_letter"] <- "C"
tmp_S15_GL[tmp_S15_GL$allele2 == 2,"a2_letter"] <- "G"
tmp_S15_GL[tmp_S15_GL$allele2 == 3,"a2_letter"] <- "T"

tmp_S15_t4 <- SUPER_15_freq$type4_freq[S15_p1_daf_filter]
tmp_S15_t1 <- SUPER_15_freq$type1_freq[S15_p1_daf_filter]

tmp_S15_t1[tmp_S15_GL$a1_letter == tmp_S15_freq$major] <- 1 - tmp_S15_t1[tmp_S15_GL$a1_letter == tmp_S15_freq$major]
tmp_S15_t4[tmp_S15_GL$a1_letter == tmp_S15_freq$major] <- 1 - tmp_S15_t4[tmp_S15_GL$a1_letter == tmp_S15_freq$major]

#All SNPs in region
tmp_S15_full_freq <- SUPER_15_0.05_freq[SUPER_15_0.05_freq$position > 11.8e6 &  SUPER_15_0.05_freq$position < 12.2e6,]

#Canadian, and other populations
pdf(file = "~/Projects/Eel/doc/SUPER_15_pop_freq.pdf", height = 8, width = 9)
par(xpd = NA, mar =c(5,4,3,10))
plot(x = SUPER_15_0.05_freq$position[SUPER_15_0.05_freq$position %in% SUPER_15_GL_clean$pos[S15_p1_daf_filter]], y = tmp_S15_freq[,"ACA"], type  = "n", ylab = "Alternative allele frequency", xlab = " Chr 15")
for(i in 9:19){
  points(x = tmp_S15_freq$position, y = tmp_S15_freq[,i], col = i - 8, pch = i-8, cex = 0.6)
}
legend(x = par("usr")[2] + 5e3, y = par("usr")[4]-0.3, col = (9:19) - 8, pch = (9:19)-8, cex = 1.2, legend = names(tmp_S15_freq[,9:19]), text.font = c(rep(1, 4), 2, rep(1, 6)))

#par(xpd = FALSE, mar =c(5,4,3,2))
plot(x = SUPER_15_0.05_freq$position[SUPER_15_0.05_freq$position %in% SUPER_15_GL_clean$pos[S15_p1_daf_filter]], y = tmp_S15_freq[,"ACA"], type  = "n", ylab = "Alternative Allele frequency", xlab = " Chr 15")
points(x = SUPER_15_GL_clean$pos[S15_p1_daf_filter], y = tmp_S15_t4, pch = 16, col = "cornflowerblue", cex = 0.4)
points(x = SUPER_15_GL_clean$pos[S15_p1_daf_filter], y = tmp_S15_t1, pch = 16, col = "darkorange", cex = 0.4)
points(x = tmp_S15_freq$position, y = 1-as.numeric(tmp_S15_GL$a2_letter == tmp_S15_freq$major), col = "black", pch = 20)
legend(x = par("usr")[1], y = 0.6, col = c("cornflowerblue", "darkorange", "black"), pch = 16, cex = 1.5, legend = c("Common (European)", "Rare (American)", "Reference equal to \"Common\""))

plot(x = tmp_S15_full_freq$position, y = tmp_S15_full_freq[,"ACA"], ylab = "Alt allele frequency", xlab = " Chr 15", main = "ACA", pch = 16, cex = 0.6)
plot(x = tmp_S15_full_freq$position, y = tmp_S15_full_freq[,"ASRI"], ylab = "Alt allele frequency", xlab = " Chr 15", main = "ASRI", pch = 16, cex = 0.6)

dev.off()



#Inset low-diveristy region
S15_rep_GR <- GRanges("SUPER_15:12123000-12140000")
writeXStringSet(filepath = "~/Projects/Eel/data/target_regions/S15_12123_12140_kb.fa", x = fAngAng1_pri[S15_rep_GR])
#Starting low-diveristy region
S15_flank_GR <- GRanges("SUPER_15:11910000-11927000")
writeXStringSet(filepath = "~/Projects/Eel/data/target_regions/S15_11910_11927_kb.fa", x = fAngAng1_pri[S15_flank_GR])

#SUPER 15 peak 2
#SUPER_15_33455012 - SUPER_15_33674630
#zgrep -m1 -n -o SUPER_15_33455012 SUPER_19_SUPER_15.beagle.gz
767893:SUPER_15_33455012
780419:SUPER_15_33674630
zcat SUPER_19_SUPER_15.beagle.gz | awk 'NR == 1 || (NR >=767893 && NR <= 780419) {print}' > SUPER_15_peak2.beagle 
zcat SUPER_19_SUPER_15.beagle.gz | awk 'NR == 1 || (NR >=747893 && NR <= 800419) {print}' > SUPER_15_peak2_ext.beagle 



#SUPER_15_peak2_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak2.beagle.gz", header = T, sep = "\t", stringsAsFactors = F)
SUPER_15_peak2_GL <- read.table(file = "~/Projects/Eel/data/geno_likelihood/SUPER_15_peak2_ext.beagle.gz", header = T, sep = "\t", stringsAsFactors = F) #Adding some flanks for context
SUPER_15_peak2_GL[,"pos"] <- as.numeric(sub("SUPER_15_", "", SUPER_15_peak2_GL[,"marker"]))
SUPER_15_peak2_GL[,"pc1.pval"] <- SUPER_15_site_sel[match(SUPER_15_peak2_GL[,"pos"], SUPER_15_site_sel[,"position"]),"pc1.pval"]

p2_MajorHomCols <- grep("Ind[0-9]+$", names(SUPER_15_peak2_GL))DE
p2_pc1_filter <- SUPER_15_peak2_GL$pc1.pval < 1e-4
p2_pc1_filter[is.na(p2_pc1_filter)] <- F

SUPER_15_p2_GL_clean <- SUPER_15_peak2_GL
for(i in 1:length(p2_MajorHomCols)){
  ind_SDs <- rowSds(as.matrix(SUPER_15_p2_GL_clean[, p2_MajorHomCols[i] + c(0,1,2)]), na.rm = T)
  SUPER_15_p2_GL_clean[ind_SDs == 0,  p2_MajorHomCols[i] + c(0,1,2)] <- NA
}

SUPER_15_p2_kinship <- region_kinship(SUPER_15_p2_GL_clean[p2_pc1_filter,], Beagle_ind_list)
SUPER_15_p2_kinship[[1]][is.na(SUPER_15_p2_kinship[[1]])] <- mean(SUPER_15_p2_kinship[[1]], na.rm = T) #AMAR10 is very problematic in this region
SUPER_15_p2_hclust <- hclust(as.dist(SUPER_15_p2_kinship[[1]]))
plot(SUPER_15_p2_hclust)
heatmap(SUPER_15_p2_kinship[[1]], scale = "none", Rowv = as.dendrogram(SUPER_15_p2_hclust), as.dendrogram(SUPER_15_p2_hclust))

#SUPER_13 1~8 Mb
#Outputting bed-file for filtering RNAseq vcfs (of males only)
SUPER_13_ext_GL$chr <- "SUPER_13"
write.table(SUPER_13_ext_GL[,c("chr", "pos")], sep = "\t", row.names = F, quote = F, col.names = F, file = "~/Projects/Eel/data/geno_likelihood/SUPER_13_block1_ext.bed")
#vcftools --gzvcf /proj/snic2020-2-19/private/eel/users/erik/rna_mapping/SRR1564733/SRR1564733_fil.vcf.gz --positions ./SUPER_13_block1_ext.bed --out ./SRR1564733 --recode-INFO-all --recode
SRR1564733 <- read.table("~/Projects/Eel/data/target_regions/SUPER_13/SRR1564733.recode.vcf", stringsAsFactors = F)
SRR1564732 <- read.table("~/Projects/Eel/data/target_regions/SUPER_13/SRR1564732.recode.vcf", stringsAsFactors = F)
vcf_files <- grep("vcf", dir("~/Projects/Eel/data/target_regions/SUPER_13", full.names = T), value = T)
vcf_names <- sub(".+(SRR15[0-9]+).recode.+", "\\1", vcf_files)
vcf_list <- list()
for(i in 1:length(vcf_files)){
  vcf_list[[vcf_names[i]]] <- read.table(vcf_files[i], stringsAsFactors = F)
}
for(i in 1:length(vcf_files)){
  vcf_list[[paste0(vcf_names[i], "_diag")]] <- vcf_list[[i]][vcf_list[[i]]$V2 %in% SUPER_13_ext_GL$pos[diag_filter],c(1:5,7,10)]
}

pos_vec <- integer()
for(i in 1:length(vcf_files)){
  pos_vec <- c(pos_vec, vcf_list[[paste0(vcf_names[i], "_diag")]][,2])
}

SUPER_13_SNP_GR <- GRanges(seqnames = "SUPER_13", ranges = IRanges(start = SUPER_13_ext_GL_clean$pos, end = SUPER_13_ext_GL_clean$pos))
SUPER_13_SNP_exon_hits <- findOverlaps(SUPER_13_SNP_GR, fAngAng_gff_rename[fAngAng_gff_rename$type %in% c("exon", "lnc_RNA")])#[fAngAng_gff_rename$type == "exon"])
exonic_pos <- SUPER_13_SNP_GR[unique(SUPER_13_SNP_exon_hits@from)]@ranges@start
diag_exonic_pos <- exonic_pos[exonic_pos %in% SUPER_13_ext_GL_clean$pos[diag_filter]]



#SUPER_1 ~81 Mb
#Approximate start of region:
zgrep -m1 -n "81190048" SUPER_1_0.05_snpEff.out.gz
2763328:SUPER_1 81190048  
zcat SUPER_1_0.05_snpEff.out.gz | awk 'NR == 1 || (NR >=2763328-1e5 && NR <= 2763328+1e5) {print}' > SUPER_1_81Mb.ann
SUPER_1_81_snpEff <- read.table("~/Projects/Eel/data/annotation/SUPER_1_81Mb.ann.gz", stringsAsFactors = F, sep = "\t")
SUPER_1_81MB_GR <- GRanges(seqnames = "SUPER_1", ranges = IRanges(start = 8.1175e7, end = 8.125e7))
load("~/Projects/Eel/data/all_maf/maf_co_0.1/SUPER_1_freq.RData")
t_daf <- abs(rowMeans(SUPER_1_freq[,midAtl_ids], na.rm = T) - rowMeans(SUPER_1_freq[,non_midAtl_ids], na.rm = T))
t_daf[is.na(t_daf)] <- 0
pos_filter <- t_daf > 0.10 & SUPER_1_freq$position > SUPER_1_81MB_GR@ranges@start & SUPER_1_freq$position < SUPER_1_81MB_GR@ranges@start + SUPER_1_81MB_GR@ranges@width & SUPER_1_freq$nInd > 300
range_filter <- SUPER_1_freq$position > SUPER_1_81MB_GR@ranges@start & SUPER_1_freq$position < SUPER_1_81MB_GR@ranges@start + SUPER_1_81MB_GR@ranges@width  & SUPER_1_freq$nInd > 300
t_snpEff <- SUPER_1_81_snpEff$V2 %in% SUPER_1_freq$position[pos_filter]
which.max(t_daf[pos_filter])
SUPER_1_freq$position[pos_filter][45]

SUPER_1_hits <- findOverlaps(SUPER_1_81MB_GR, fAngAng_gff_rename)
SUPER_1_81_gff <- fAngAng_gff_rename[SUPER_1_hits@to]
SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]
SUPER_1_81_miss_pos <- SUPER_1_81_snpEff[t_snpEff,2][grep("missense", SUPER_1_81_snpEff$V8[t_snpEff])]
SUPER_1_81_miss_DAF <- t_daf[match(SUPER_1_81_miss_pos, SUPER_1_freq$position)]

pdf(file = "~/Projects/Eel/doc/SUPER_1_81MB_annot.pdf", width = 10)
plot(x= SUPER_1_freq$position[range_filter], y = t_daf[range_filter], pch = 20, cex = 0.5, col = "grey70", xlab = "Position", ylab = "DAF", ylim = c(0, 0.22), xlim = c(8.1180e7,8.1205e7)) 
segments(x0 = SUPER_1_81_miss_pos, y0 = SUPER_1_81_miss_DAF, y1 = 0.2, lty = "dashed", col = "grey50") #Missense location
#abline(v = 81196156) #Peak location
points(x= SUPER_1_freq$position[pos_filter], y = t_daf[pos_filter], pch = 20, col = "cornflowerblue", cex = 0.8)
segments(x0 = SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]@ranges@start[2], x1 = (SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]@ranges@start+SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]@ranges@width)[2], y0 = 0.20, col = "firebrick", lwd = 3)
x0_vec <- SUPER_1_81_gff[SUPER_1_81_gff$type == "exon"]@ranges@start
x1_vec <- x0_vec + SUPER_1_81_gff[SUPER_1_81_gff$type == "exon"]@ranges@width
gene_mid_vec <- SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]@ranges@start+ (0.5*SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]@ranges@width)
rect(xleft = x0_vec, xright = x1_vec, ybottom = 0.19, ytop = 0.21, col = "firebrick", border = NA)
text(x = gene_mid_vec[2], y = 0.22, labels =  SUPER_1_81_gff[SUPER_1_81_gff$type == "gene"]$Name[2], cex =1.5)
points(x = SUPER_1_81_miss_pos , y = SUPER_1_81_miss_DAF, col = "red", cex = 1, pch = 16)
dev.off()

#SUPER_16, 21.66 MB
SUPER_16_21.6_MB_GR <- GRanges(seqnames = "SUPER_16", ranges = IRanges(start = 2.16e7, end = 2.19e7))
load("~/Projects/Eel/data/all_maf/maf_co_0.05/SUPER_16_0.05_freq.RData")
S16_21.6_sites <- SUPER_16_0.05_freq$position > 2.16e7 & SUPER_16_0.05_freq$position < 2.19e7 & SUPER_16_0.05_freq$nInd > 300
SUPER_16_22MB_freq <- SUPER_16_0.05_freq[S16_21.6_sites,-(1:8)]
rownames(SUPER_16_22MB_freq) <- SUPER_16_0.05_freq[S16_21.6_sites,"position"]
#high_daf <- abs(rowMeans(SUPER_16_22MB_freq[,c("AAF", "AI", "AES")]) - rowMeans(SUPER_16_22MB_freq[,c("ASRI", "SMS", "AST")])) > 0.2
SUPER_16_22MB_DAF <- abs(rowMeans(SUPER_16_22MB_freq[,c("AAF", "AES")]) - rowMeans(SUPER_16_22MB_freq[,c("ASRI", "SMS", "AST")]))
S16_high_daf <- SUPER_16_22MB_DAF > 0.2

SUPER_16_hits <- findOverlaps(SUPER_16_21.6_MB_GR, fAngAng_gff_rename)
SUPER_16_21.6_gff <- fAngAng_gff_rename[SUPER_16_hits@to]
SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]

SUPER_16_snpEff <- read.table("~/Projects/Eel/data/annotation/SUPER_16_0.05_snpEff.out.gz", stringsAsFactors = F, sep = "\t")
SUPER_16_22MB_snpEff <- SUPER_16_snpEff[SUPER_16_snpEff$V2 %in% rownames(SUPER_16_22MB_freq),]
SUPER_16_22MB_miss <- grep("missense", SUPER_16_22MB_snpEff$V8)

pdf(file = "~/Projects/Eel/doc/SUPER_16_22MB_annot.pdf")
heatmap(t(as.matrix(SUPER_16_22MB_freq[S16_high_daf,])), scale = "none", Colv = NA)
plot(y = SUPER_16_22MB_DAF, x = rownames(SUPER_16_22MB_freq), pch = 16, xlab = "Position on SUPER_16", main = "AES + AAF v SWE", ylab = "DAF", ylim = c(0,0.45), xlim = c(21.6e6, 21.75e6) ) 
segments(x0 = SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]@ranges@start, x1 = (SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]@ranges@start+SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]@ranges@width), y0 = 0.40, col = "firebrick", lwd = 3)
x0_vec <- SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "exon"]@ranges@start
x1_vec <- x0_vec + SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "exon"]@ranges@width
gene_mid_vec <- SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]@ranges@start + (0.5*SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]@ranges@width)
rect(xleft = x0_vec, xright = x1_vec, ybottom = 0.39, ytop = 0.41, col = "firebrick", border = NA)
points(x = SUPER_16_22MB_snpEff$V2[SUPER_16_22MB_miss], y = SUPER_16_22MB_DAF[SUPER_16_22MB_miss], col = "red", cex = 1.2, pch = 16)
text(x = gene_mid_vec, y = c(0.42, 0.38), labels =  SUPER_16_21.6_gff[SUPER_16_21.6_gff$type == "gene"]$Name, cex =0.8)
dev.off()

#Support functions
sliding_deviation_v2 <- function(reg_GL, clust_vec, ref_freq_vec, pdf_file = NULL, png_file = NULL, col_vec = NULL, smooth_int = 100, target_gff = NULL){
  MajorHomCols <- grep("Ind[0-9]+$", names(reg_GL))
  lwd_val <- 0.3 
  if(!is.null(png_file)) lwd_val <- 1
  
  if(is.null(col_vec[1])){
    col_vec <- col2rgb(c("red", "darkorchid", "grey50", "blue"),alpha = T)
    col_vec["alpha",] <- 80
  }
  
  if(!is.null(pdf_file)) pdf(pdf_file, width = 10)
  if(!is.null(png_file)) png(png_file, width = 2000, height = 1500)
  plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = paste0("Allele count deviaiton from reference group mean(",smooth_int," SNP rolling avg)"), main = "", xlab = "Position")
  for(i in 1:length(MajorHomCols)){ #length(MajorHomCols)
    kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - 2*ref_freq_vec)
    #Adjusting for the group avearge at each marker, and also offsetting according to the mean
    if(sum(!is.na(kinship_dev_vec)) > smooth_int){
      lines(x = reg_GL[,"pos"][!is.na(kinship_dev_vec)], y = filter(kinship_dev_vec[!is.na(kinship_dev_vec)], rep(1/smooth_int, smooth_int)), lty = (i %% 3) +1, lwd = lwd_val, col = rgb(t(col_vec[,clust_vec[i]]), maxColorValue = 255))
    }
  }
  if(!is.null(target_gff[1])){
    segments(x0 = target_gff[target_gff$type == "gene"]@ranges@start, x1 = target_gff[target_gff$type == "gene"]@ranges@start+target_gff[target_gff$type == "gene"]@ranges@width, y0 = 1.6, col = "firebrick", lwd = 3)
    
    x0_vec <- target_gff[target_gff$type == "exon"]@ranges@start
    x1_vec <- x0_vec + target_gff[target_gff$type == "exon"]@ranges@width
    rect(xleft = x0_vec, xright = x1_vec, ybottom = 1.57, ytop = 1.63, col = "firebrick", border = NA)
    gene_mid_vec <- target_gff[target_gff$type == "gene"]@ranges@start+ (0.5*target_gff[target_gff$type == "gene"]@ranges@width)
    text(x = gene_mid_vec, y = c(1.65, 1.55), labels = target_gff[target_gff$type == "gene"]$Name, cex = 1.5)
  }
  
  if(!is.null(png_file)) dev.off()
  #col_vec_hl <- col_vec
  #col_vec_hl["alpha",] <- 150
  
  for(j in unique(clust_vec)){
    if(!is.null(png_file)) png(paste0(sub(".png", "", png_file), "_type_", j, ".png"), width = 2000, height = 1500)
    plot(x= reg_GL[,"pos"], y = reg_GL[,MajorHomCols[1]], type = "n", ylim = c(0,2), ylab = paste0("Allele count deviaiton from reference group mean(",smooth_int," SNP rolling avg)"), main = "", xlab = "Position")
    for(i in which(clust_vec != j)){ #length(MajorHomCols)
      kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - 2*ref_freq_vec)
      #Adjusting for the group avearge at each marker, and also offsetting according to the mean
      if(sum(!is.na(kinship_dev_vec)) > smooth_int){
        lines(x = reg_GL[,"pos"][!is.na(kinship_dev_vec)], y = filter(kinship_dev_vec[!is.na(kinship_dev_vec)], rep(1/smooth_int, smooth_int)), lty = (i %% 3) +1, lwd = lwd_val, col = "grey70") 
      }
    }
    for(i in which(clust_vec == j)){ #length(MajorHomCols)
      kinship_dev_vec <- abs((2*reg_GL[,MajorHomCols[i]] + reg_GL[,MajorHomCols[i]+1]) - 2*ref_freq_vec)
      #Adjusting for the group avearge at each marker, and also offsetting according to the mean
      if(sum(!is.na(kinship_dev_vec)) > smooth_int){
        lines(x = reg_GL[,"pos"][!is.na(kinship_dev_vec)], y = filter(kinship_dev_vec[!is.na(kinship_dev_vec)], rep(1/smooth_int, smooth_int)), lty = (i %% 3) +1, lwd = lwd_val, col = rgb(t(col_vec[,clust_vec[i]]), maxColorValue = 255)) 
      }
    }
    if(!is.null(png_file)) dev.off()
  }
  if(!is.null(pdf_file)) dev.off()
}
