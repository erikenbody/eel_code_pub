#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 8 #needs full node for array 1-6
#SBATCH -t 3-00:00:00
#SBATCH -M snowy
#SBATCH -J angsd
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out
module load samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd_2Sep2020/angsd
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd

WORK_D=/proj/snic2020-2-19/private/eel/users/erik/thetas/super_15_taj_d_calcs

REFGENOME=/proj/uppoff2019007/private/eel/reference/fAngAng1.pri.cur.20200204.fasta
cd $WORK_D

INTERVAL=SUPER_15

POP_BAMLIST=super_15_1.txt

OUTNAME=group1_${INTERVAL}
$ANGSD_PATH/angsd -bam $POP_BAMLIST -r $INTERVAL -doSaf 1 -anc $REFGENOME -GL 2 -P 8 -out $OUTNAME -doCounts 1 -setMinDepth 5 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 1 -minQ 20 -remove_bads 1 -uniqueOnly 1 -dumpCounts 2

$ANGSD_PATH/misc/realSFS ${OUTNAME}.saf.idx -P 8 -fold 1 > ${OUTNAME}.sfs
$ANGSD_PATH/misc/realSFS saf2theta  ${OUTNAME}.saf.idx -sfs ${OUTNAME}.sfs -outname ${OUTNAME}
$ANGSD_PATH/misc/thetaStat do_stat ${OUTNAME}.thetas.idx -win 10000 -step 10000  -outnames ${OUTNAME}.theta.thetasWindow.gz
$ANGSD_PATH/misc/thetaStat do_stat ${OUTNAME}.thetas.idx -win 1000 -step 1000  -outnames ${OUTNAME}.theta.1kb.thetasWindow.gz

##
POP_BAMLIST=super_15_4.txt

OUTNAME=group4_${INTERVAL}
$ANGSD_PATH/angsd -bam $POP_BAMLIST -r $INTERVAL -doSaf 1 -anc $REFGENOME -GL 2 -P 8 -out $OUTNAME -doCounts 1 -setMinDepth 5 -setMaxDepth 1000 -setMinDepthInd 0.25 -minMapQ 1 -minQ 20 -remove_bads 1 -uniqueOnly 1 -dumpCounts 2

$ANGSD_PATH/misc/realSFS ${OUTNAME}.saf.idx -P 8 -fold 1 > ${OUTNAME}.sfs
$ANGSD_PATH/misc/realSFS saf2theta  ${OUTNAME}.saf.idx -sfs ${OUTNAME}.sfs -outname ${OUTNAME}
$ANGSD_PATH/misc/thetaStat do_stat ${OUTNAME}.thetas.idx -win 10000 -step 10000  -outnames ${OUTNAME}.theta.thetasWindow.gz
$ANGSD_PATH/misc/thetaStat do_stat ${OUTNAME}.thetas.idx -win 1000 -step 1000  -outnames ${OUTNAME}.theta.1kb.thetasWindow.gz
