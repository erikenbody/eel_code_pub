#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se




#Eel depth data
#SUPER_13 region of interest (PCA-based)
#Find the start of SUPER_13 in the pos-file



zgrep -n -o -m1 SUPER_13 /crex/proj/uppstore2017191/private/users/erik/eels/thetas/no_maf_cutoff_files/euro.pos.gz
266726808:SUPER_13
#extract every 1000th position for the next 1e7 lines
zcat /crex/proj/uppstore2017191/private/users/erik/eels/thetas/no_maf_cutoff_files/euro.pos.gz | awk 'NR == 1 || (NR >=266726808 && NR <= 266726808+1e7 && NR % 1000 == 1) {print}' > SUPER_13_thinned_out_depth.pos 
zcat /crex/proj/uppstore2017191/private/users/erik/eels/thetas/no_maf_cutoff_files/euro.counts.gz | awk 'NR == 1 || (NR >=266726808 && NR <= 266726808+1e7 && NR % 1000 == 1) {print}' > SUPER_13_thinned_out_depth.ind

SUPER_13_dp <- read.table("~/Projects/Eel/data/depth/SUPER_13_thinned_out_depth.pos", stringsAsFactors = F, header = T, sep = "\t")
SUPER_13_dp_ind <- read.table("~/Projects/Eel/data/depth/SUPER_13_thinned_out_depth.ind", stringsAsFactors = F, header = T, sep = "\t")
SUPER_13_dp_ind <- cbind(SUPER_13_dp, SUPER_13_dp_ind[,-dim(SUPER_13_dp_ind)[2]]) #There is an empty column at the end of the file 

class2_cols <- which(Beagle_ind_list$super_13_clust_v2 == 2) + 3 #Offsetting the initial columns
class1_cols <- which(Beagle_ind_list$super_13_clust_v2 == 1) + 3 #Offsetting the initial columns

SUPER_13_GR <- GRanges(seqnames = (fAngAng_gff[fAngAng_gff$type == "region"][13])@seqnames, ranges = IRanges(start = 1, end = 1e7))

SUPER_13_hits <- findOverlaps(SUPER_13_GR, fAngAng_gff)
SUPER_13_gff <- fAngAng_gff[SUPER_13_hits@to]

pdf(file = "~/Projects/Eel/doc/SUPER_13_depth.pdf", width = 10, height = 7)
plot(x = SUPER_13_dp_ind$pos, y = filter(rowSums(SUPER_13_dp_ind[,class1_cols])/length(class1_cols), rep(1/20,20)), pch = 16, col = "grey70", ylim = c(0,3), cex  = 0.3)
points(x = SUPER_13_dp_ind$pos, y = filter(rowSums(SUPER_13_dp_ind[,class2_cols])/length(class2_cols), rep(1/20,20)), pch = 16, col = "grey30", ylim = c(0,3), cex = 0.3)

points(x = SUPER_13_ext_GL_clean$pos[!diag_filter], y = abs(SUPER_13_ext_freq$type2_freq[!diag_filter]-SUPER_13_ext_freq$type1_freq[!diag_filter]), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3)
points(x = SUPER_13_ext_GL_clean$pos[diag_filter], y = abs(SUPER_13_ext_freq$type2_freq[diag_filter]-SUPER_13_ext_freq$type1_freq[diag_filter]), pch = 16, col = "darkorchid", cex = 0.3)

segments(x0 = SUPER_13_gff[SUPER_13_gff$type == "gene"]@ranges@start, x1 = SUPER_13_gff[SUPER_13_gff$type == "gene"]@ranges@start+SUPER_13_gff[SUPER_13_gff$type == "gene"]@ranges@width, y0 = 2.5, col = "firebrick", lwd = 3)
dev.off()

#SUPER_15
#zgrep -n -o -m1 SUPER_15 /crex/proj/uppstore2017191/private/users/erik/eels/thetas/no_maf_cutoff_files/euro.pos.gz
#MIssing form depth file "unexpected EOF"

#Call rate depth estiamtion (p_0 = e^(-cov) => cov = -ln(p_0))
pdf(file = "~/Projects/Eel/doc/SUPER_15_depth.pdf", width = 10)
plot(x = SUPER_15_GL_clean$pos[pc1_filter], y = -log(rowSums(is.na(SUPER_15_GL_clean[pc1_filter,MajorHomCols]))/length(MajorHomCols)), pch = 16, col = "darkorchid", cex = 0.3, ylim = c(0,5), main = "SUPER_15, peak 1", xlab = "Position", ylab = "Estimated depth")
points(x = SUPER_15_GL_clean$pos[!pc1_filter], y = -log(rowSums(is.na(SUPER_15_GL_clean[!pc1_filter,MajorHomCols]))/length(MajorHomCols)), pch = 16, col = "grey50", cex = 0.3)
legend(x = "topleft", legend = c("Background markers","Diagnostic markers"), col  = c("grey50", "darkorchid"), pch = 16)

plot(x = SUPER_15_p2_GL_clean$pos[!p2_pc1_filter], y = -log(rowSums(is.na(SUPER_15_p2_GL_clean[!p2_pc1_filter,p2_MajorHomCols]))/length(p2_MajorHomCols)), pch = 16, col = "grey50", ylim = c(0, 5), cex = 0.3, main = "SUPER_15, peak 2", xlab = "Position", ylab = "Estimated depth")
points(x = SUPER_15_p2_GL_clean$pos[p2_pc1_filter], y = -log(rowSums(is.na(SUPER_15_p2_GL_clean[p2_pc1_filter,p2_MajorHomCols]))/length(p2_MajorHomCols)), pch = 16, col = "darkorchid", cex = 0.3)
legend(x = "topleft", legend = c("Background markers","Diagnostic markers"), col  = c("grey50", "darkorchid"), pch = 16)
dev.off()

