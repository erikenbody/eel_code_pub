#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
data_dir=/proj/snic2020-2-19/private/eel/users/erik/rna_mapping
bam_in=$1

# Update with the location of the reference data files
fasta=$2

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/Uppsala_University_eval.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/sentieon-genomics-201911

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
sample=$3
group=$4
platform="ILLUMINA"

# Other settings
nt=10 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$data_dir/$sample #ede edited this to be sample
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i $bam_in --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i $bam_in --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option ${sample}.deduped.bam
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -r $fasta -i ${sample}.deduped.bam \
  --algo RNASplitReadsAtJunction --reassign_mapq 255:60 ${sample}.split.bam
# ******************************************
# 6. HC Variant caller
# ******************************************

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -r $fasta -i ${sample}.split.bam \
  --algo Haplotyper --trim_soft_clip  \
  --call_conf 20 --emit_conf 20 ${sample}.vcf
