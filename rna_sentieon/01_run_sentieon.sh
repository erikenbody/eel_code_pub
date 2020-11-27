#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 10
#SBATCH -M snowy
#SBATCH -t 5-00:00:00
#SBATCH -J eel_sentieon
#SBATCH -e eel_sentieon_%A_%a.err            # File to which STDERR will be written
#SBATCH -o eel_sentieon_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml python/2.7.15 star/2.7.2b samtools/1.10 sratools/2.9.6-1

TOPDIR=/proj/snic2020-2-19/private/eel/users/erik/rna_mapping
REFERENCE=/proj/snic2020-2-19/private/eel/users/mats/reference/fAngAng1.pri.cur.20200204.fasta
READS=/proj/snic2020-2-19/private/eel/reads/RNA_reads/SRA
BAMADDRG=/crex/proj/snic2020-2-19/private/eel/users/erik/rna_mapping/bamaddrg
cd $TOPDIR

NAME=`cat $TOPDIR/eel_sra_list.txt| awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPBASE=$(basename $NAME)
SAMPLE=${SAMPBASE%.sra*}

##split output sra
cd $READS
fastq-dump -I --gzip --split-files $NAME
cd $TOPDIR

##set up naming
FILENAME=$READS/${SAMPLE}_1.fastq.gz

SAMPBASE=$(basename $FILENAME)
SAMPLE=${SAMPBASE%_1*}
SAMPLE_DIR=$(dirname $FILENAME)
R1=$SAMPLE_DIR/${SAMPLE}*_1*fastq.gz
R2=$SAMPLE_DIR/${SAMPLE}*_2*fastq.gz
mkdir -p results

##run star aligner following gatk4 best practices
STAR --genomeDir /crex/proj/snic2020-2-19/private/eel/users/erik/rna_mapping/eel_genome_idx \
--runThreadN 10 \
--readFilesIn $R1 $R2 \
--readFilesCommand "gunzip -c" \
--outFileNamePrefix ./results/${SAMPLE}_ \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate

##add a read group header, could this this by star if you re write this
$BAMADDRG/bamaddrg -b ./results/${SAMPLE}_Aligned.sortedByCoord.out.bam -s $SAMPLE -r sra.eels > ./results/${SAMPLE}_Aligned.sortedByCoord.out_rh.bam

BAM=`realpath ./results/${SAMPLE}_Aligned.sortedByCoord.out_rh.bam`
samtools index $BAM

##run sentieon haplotype caller
chmod +x ~/ec/rna_sentieon/eel_sentieon.sh
sh ~/ec/rna_sentieon/eel_sentieon.sh $BAM $REFERENCE $SAMPLE $SAMPBASE

##run gatk variant filter best practicese
ml GATK/4.1.4.1

gatk \
     VariantFiltration \
   --R $REFERENCE \
   --V ${SAMPLE}/${SAMPLE}.vcf \
   --window 35 \
   --cluster 3 \
   --filter-name "FS" \
   --filter-expression "FS > 30.0" \
   --filter-name "QD" \
   --filter-expression  "QD < 2.0" \
   -O ${SAMPLE}/${SAMPLE}_fil.vcf.gz

gzip ${SAMPLE}/${SAMPLE}.vcf
