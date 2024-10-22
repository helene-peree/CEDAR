#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=multiQC2
#SBATCH --output=../slurm_logs/5_multiQC_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

set -e
set -u
set -o pipefail

date
echo -e "\n"
echo "node used: $SLURM_JOB_NODELIST"
echo "memory per CPU: $SLURM_MEM_PER_CPU"
echo -e "\n"

module load EasyBuild
#module load MultiQC/0.9-foss-2016b-Python-2.7.12
module load MultiQC/1.14-foss-2022b

Plate=$1
Software=$2 #STAR or WASP

WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene"
InOutputDir=${WorkingDir}/RNAseq/aln/${Plate}/${Software}_QC

srun multiqc ${InOutputDir}/*RNAseq_metrics.out --force --filename ${InOutputDir}/0_multiqc_report_RNAseq_metrics.html

echo -e "\n"
echo "QC analysis with multiQC completed"
echo -e "\n"
date
