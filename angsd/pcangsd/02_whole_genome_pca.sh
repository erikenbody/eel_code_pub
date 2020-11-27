#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p node
#SBATCH -t 5-00:00:00
#SBATCH -J angsd_pca
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_pca_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_pca_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd_2Sep2020/angsd
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd

WORK_D=/home/eenbody/snic2020-2-19/private/eel/users/erik/pcangsd_overflow
REFGENOME=/proj/uppoff2019007/private/eel/reference/fAngAng1.pri.cur.20200204.fasta
BAMLIST=/crex/proj/uppstore2017191/private/users/erik/eels/eel_bams.txt
cd $WORK_D
SITES_FILE=/crex/proj/snic2020-2-19/private/eel/users/mats/site_files/all_eel_scaffs_maf_0.1_sites.txt
#
$ANGSD_PATH/angsd -ref $REFGENOME -sites $SITES_FILE \
          -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
          -GL 2 -out eels_all_V2 -nThreads 20 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -bam $BAMLIST

pyenv global 3.7.0

python $PCANGSD_PATH/pcangsd.py -beagle eels_all_V2.beagle.gz -o eels_all_V2 -admix -threads 20
