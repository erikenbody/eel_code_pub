#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Annotation
require(rtracklayer)
fAngAng_gff <- rtracklayer::import("~/Projects/Eel/data/annotation/GCF_013347855.1_fAngAng1.pri_genomic.gff.gz", format = "GFF")
require(Biostrings)
#Reference genome
fAngAng1_pri <- readDNAStringSet("~/Projects/Eel/data/reference/fAngAng1.pri.cur.20200204.fasta.gz")

#Annotaion are ordered by size, matching the SUPER_number
fAngAng_gff[fAngAng_gff$type == "region"][13]
fAngAng1_pri["SUPER_13"]

#Attempting to rename the chromosomes in the annotation
fAngAng_gff_rename <- fAngAng_gff

ref_seqname_df <- as.data.frame(fAngAng1_pri@ranges)[order(fAngAng1_pri@ranges@width, decreasing = T),]
ann_seqname_df <- as.data.frame(fAngAng_gff[fAngAng_gff$type == "region"])[-54,1:5] #removing the mitochodria entry
ref_seqname_df <- ref_seqname_df[match(ann_seqname_df$width, ref_seqname_df$width),]
joint_seqname_df <- cbind(ref_seqname_df, ann_seqname_df)

seqinfo(fAngAng_gff_rename, 1:54) <- Seqinfo(seqnames = c("mt", joint_seqname_df$names))
fAngAng_gff_rename[fAngAng_gff_rename$type == "region"]@seqnames[-54] == joint_seqname_df$names
fAngAng_gff_rename[fAngAng_gff_rename$type == "region"]@ranges@width[-54] == joint_seqname_df[,4] 
#Checks out.

rtracklayer::export(fAngAng_gff_rename, con = "~/Projects/Eel/data/annotation/fAngAng1.gff", format = "gff3")
awk -v OFS='\t' 'NR == 1 {printf "%s", "#"}; {print $1, $2, $1"_"$2, $5, $4}' snpEff_test2.txt #Reordering the columns to match snpEFF expectations
snpEff ann -v -dataDir ./ -c ./snpEff.config fAngAng1 SUPER_1_maf_0.05_reformat.txt > SUPER_1_0.05_snpEff.out #The actual annotation step

SUPER_15_snpEff <- read.table("~/Projects/Eel/data/annotation/SUPER_15_high_DAF_snpEff.out", stringsAsFactors = F, sep = "\t")
SUPER_15_snpEff_high_DAF <- SUPER_15_snpEff[SUPER_15_snpEff$V2 %in% SUPER_15_GL_clean$pos[pc1_filter],]

SUPER_15_GR <- GRanges(seqnames = (fAngAng_gff[fAngAng_gff$type == "region"][15])@seqnames, ranges = IRanges(start = 11.8e6, end = 12.2e6))
SUPER_15_hits <- findOverlaps(SUPER_15_GR, fAngAng_gff)
SUPER_15_gff <- fAngAng_gff[SUPER_15_hits@to]
SUPER_15_gff[SUPER_15_gff$type == "gene"]

pdf(file = "~/Projects/Eel/doc/SUPER_15_clean_anont.pdf", width = 10)
plot(x = SUPER_15_GL_clean$pos[!pc1_filter], y = abs(SUPER_15_freq$type4_freq[!pc1_filter]-SUPER_15_freq$type1_freq[!pc1_filter]), ylim = c(0,1), pch = 16, col = "grey50", cex = 0.3, ylab = "DAF", xlab = "Position")
points(x = SUPER_15_GL_clean$pos[pc1_filter], y = abs(SUPER_15_freq$type4_freq[pc1_filter]-SUPER_15_freq$type1_freq[pc1_filter]), pch = 16, col = "darkorchid", cex = 0.3)
abline(v = SUPER_15_snpEff_high_DAF[grep("missense", SUPER_15_snpEff_high_DAF$V8),2])

segments(x0 = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start, x1 = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start+SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@width, y0 = 0.9, col = "firebrick", lwd = 3)
abline(v = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start,  col = "firebrick1", lty = "dashed", lwd = 0.5)
abline(v = SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@start + SUPER_15_gff[SUPER_15_gff$type == "gene"]@ranges@width, col = "firebrick1", lty = "dashed", lwd = 0.5)

x0_vec <- SUPER_15_gff[SUPER_15_gff$type == "exon"]@ranges@start
x1_vec <- x0_vec + SUPER_15_gff[SUPER_15_gff$type == "exon"]@ranges@width
rect(xleft = x0_vec, xright = x1_vec, ybottom = 0.88, ytop = 0.92, col = "firebrick", border = NA)
dev.off()



