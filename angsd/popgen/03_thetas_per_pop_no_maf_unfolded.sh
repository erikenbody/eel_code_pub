#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 8
#SBATCH -M snowy
#SBATCH -t 3-10:00:00
#SBATCH -J angsd
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd_2Sep2020/angsd
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd

WORK_D=/proj/snic2020-2-19/private/eel/users/erik/thetas/per_pop/unfolded
REFGENOME=/crex/proj/snic2020-2-19/private/eel/users/mats/reference/fAngAng1.pri.cur.20200204.fasta
ANC=/crex/proj/snic2020-2-19/private/eel/users/erik/thetas/ancestral_state/ACA.fa.gz

cd $WORK_D

POP_BAMLIST=`ls -1 *_bamlists.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

OUTNAME=${POP_BAMLIST/_bamlists.txt/}

$ANGSD_PATH/angsd -bam $POP_BAMLIST -doSaf 1 -anc $ANC -GL 2 -P 8 -out $OUTNAME -doCounts 1 -setMinDepth 15 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 1 -minQ 20 -remove_bads 1 -uniqueOnly 1 -dumpCounts 2 -doMajorMinor 5 -doMaf 2

##mem problem
##$ANGSD_PATH/misc/realSFS ${OUTNAME}.saf.idx -P 8 > ${OUTNAME}.sfs

while read SCAFF
do
  echo $SCAFF
  $ANGSD_PATH/misc/realSFS ${OUTNAME}.saf.idx -P 8 -r $SCAFF > ${SCAFF}_${OUTNAME}.sfs
done < scaffs.list

cat *_${OUTNAME}.sfs > ${OUTNAME}_cat.sfs

Rscript --vanilla ~/ec/angsd/popgen/script_sum_sfs.R ${OUTNAME}_cat.sfs

#https://github.com/ANGSD/angsd/issues/312
$ANGSD_PATH/misc/realSFS saf2theta ${OUTNAME}.saf.idx -sfs ${OUTNAME}_cat_rsum.sfs -outname ${OUTNAME}
$ANGSD_PATH/misc/thetaStat do_stat ${OUTNAME}.thetas.idx -win 5000 -step 5000  -outnames ${OUTNAME}.theta.thetasWindow.gz
