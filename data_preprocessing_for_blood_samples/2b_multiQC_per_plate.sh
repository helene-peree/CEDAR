#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=multiQC1
#SBATCH --output=../slurm_logs/2_multiQC_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

set -e
set -u
set -o pipefail

date
echo -e "\n"

module load EasyBuild
#module load MultiQC/0.9-foss-2016b-Python-2.7.12
module load MultiQC/1.14-foss-2022b

Plate=$1

WorkingDir='/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene'
InputDir=${WorkingDir}/RNAseq/data/fastQC/${Plate}
OutputDir=${WorkingDir}/RNAseq/data/multiQC_per_plate/${Plate}
mkdir -p ${OutputDir}

srun multiqc ${InputDir} --force --filename ${OutputDir}/multiqc_report.html

echo -e "\n"
echo "QC analysis with multiQC completed"
echo -e "\n"
date
