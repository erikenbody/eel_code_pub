#! /bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 1
#SBATCH -t 1-18:00:00
#SBATCH -J goleft
#SBATCH --mail-user erik.enbody@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -e goleft_%A_%a.err            # File to which STDERR will be written
#SBATCH -o goleft_%A_%a.out

#CAN BE RUN INTERACTIVELY

ml conda/latest
source conda_init.sh
CONDA_ENVS_PATH=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools
#conda install -c conda-forge pixy
#conda create --name pixy -c conda-forge pixy
conda activate goleft


WORK_D=/proj/snic2020-2-19/private/eel/users/erik/goleft_indexcov
cd $WORK_D

BAMPATH=/crex/proj/uppstore2017191/private/users/mats/Eel/bams

goleft indexcov --directory . $BAMPATH/*.bam
