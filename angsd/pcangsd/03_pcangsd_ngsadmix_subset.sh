#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 10
#SBATCH -t 00-20:00:00
#SBATCH -J subset_run
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e subset_run_%A_%a.err            # File to which STDERR will be written
#SBATCH -o subset_run_%A_%a.out

module load bioinfo-tools samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd_2Sep2020/angsd/
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd
REFGENOME=/proj/uppoff2019007/private/eel/reference/fAngAng1.pri.cur.20200204.fasta
BAMLIST=/crex/proj/uppstore2017191/private/users/erik/eels/eel_bams.txt

WORK_D=/home/eenbody/snic2020-2-19/private/eel/users/erik/pcangsd_overflow/subset_whole_genome
cd $WORK_D

#subset to 2 mil SNPs
FILE=/crex/proj/snic2020-2-19/private/eel/users/mats/site_files/all_eel_scaffs_maf_0.1_sites.txt
cat $FILE | awk 'NR % 10 == 0' > all_eel_scaffs_maf_0.1_sites_2mil.txt
#some will be filtered so I overshot a bit
SITES_FILE=all_eel_scaffs_maf_0.1_sites_2mil.txt
$ANGSD_PATH/angsd sites index all_eel_scaffs_maf_0.1_sites_2mil.txt

$ANGSD_PATH/angsd -ref $REFGENOME -sites $SITES_FILE \
          -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
          -GL 2 -out eels_all_subset -nThreads 16 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -bam $BAMLIST

pyenv global 3.7.0

python $PCANGSD_PATH/pcangsd.py -beagle eels_all_subset.beagle.gz -o eels_all_subset -admix -threads 16

$ANGSD_PATH/misc/NGSadmix -likes eels_all_subset.beagle.gz -K 2 -outfiles ngs_admix_k2_super1 -P 16
$ANGSD_PATH/misc/NGSadmix -likes eels_all_subset.beagle.gz -K 3 -outfiles ngs_admix_k3_super1 -P 16

###
without canada
zcat eels_all_subset.beagle.gz | cut -f 1-147,295-1485 | gzip > eels_all_subset_no_canada.beagle.gz
python $PCANGSD_PATH/pcangsd.py -beagle eels_all_subset_no_canada.beagle.gz -o eels_all_subset_no_cananda -threads 5

$ANGSD_PATH/misc/NGSadmix -likes eels_all_subset_no_canada.beagle.gz -K 2 -outfiles ngs_admix_k2_no_canada -P 5

$ANGSD_PATH/misc/NGSadmix -likes eels_all_subset_no_canada.beagle.gz -K 1 -outfiles ngs_admix_k1_no_canada -P 10
