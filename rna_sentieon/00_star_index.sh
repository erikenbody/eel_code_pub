#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 1
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J staridx
#SBATCH -e staridx%A_%a.err            # File to which STDERR will be written
#SBATCH -o staridx%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml  star/2.7.2b
module load bwa/0.7.17 picard/2.20.4

TOPDIR=/proj/snic2020-2-19/private/eel/users/erik/rna_mapping
cd $TOPDIR
mkdir -p eel_genome_idx
REFERENCE=/proj/snic2020-2-19/private/eel/users/mats/reference/fAngAng1.pri.cur.20200204.fasta
REF_DIR=$(dirname "${REFERENCE}")
REF=$(basename "${REFERENCE}")
REFBASE=${REF/.fasta/}

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir eel_genome_idx \
--genomeFastaFiles $REFERENCE


java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary \
      R=$REFERENCE \
      O=$REF_DIR/${REFBASE}.dict
