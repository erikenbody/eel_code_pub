#! /bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 10
#SBATCH -M snowy
#SBATCH -t 3-00:00:00
#SBATCH -J angsd
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e angsd_%A_%a.err            # File to which STDERR will be written
#SBATCH -o angsd_%A_%a.out

module load samtools/1.10

ANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/angsd_2Sep2020/angsd
PCANGSD_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/pcangsd

cd /crex/proj/snic2020-2-19/private/eel/users/erik/thetas/ancestral_state
POP_BAMLIST=/crex/proj/snic2020-2-19/private/eel/users/erik/thetas/per_pop/ACA_bamlists.txt
OUTNAME=ACA
$ANGSD_PATH/angsd -bam $POP_BAMLIST -dofasta 2 -P 10 -out $OUTNAME -setMinDepth 50 -doCounts 1 -setMaxDepth 1000 -minMapQ 1 -minQ 20 -remove_bads 1 -uniqueOnly 1
