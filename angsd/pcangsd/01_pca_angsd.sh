#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 10
#SBATCH -t 1-18:00:00
#SBATCH -J angsd_pca
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_pca_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_pca_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd

WORK_D=/crex/proj/uppstore2017191/private/users/erik/eels
REFGENOME=/proj/uppoff2019007/private/eel/reference/fAngAng1.pri.cur.20200204.fasta
EEL_BAMS=/crex/proj/uppstore2017191/private/users/mats/Eel/bams
INTERVAL_FILE=/crex/proj/uppstore2017191/private/users/erik/eels/superscaff.list
INTERVAL=`cat $INTERVAL_FILE | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

cd $WORK_D
realpath $EEL_BAMS/*bam > eel_bams.txt

BAMLIST=eel_bams_no_canada.txt

SITES_FILE=/crex/proj/uppstore2017191/private/users/mats/Eel/angsd_maf_wd/site_files/${INTERVAL}_maf_0.1_sites.txt

$ANGSD_PATH/angsd -ref $REFGENOME -r $INTERVAL -sites $SITES_FILE \
          -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 \
          -GL 2 -out ${INTERVAL} -nThreads 20 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -bam $BAMLIST

pyenv global 3.7.0

python $PCANGSD_PATH/pcangsd.py -beagle ${INTERVAL}.beagle.gz -o eels_${INTERVAL} -selection -sites_save -indf_save -admix -threads 10
