#!/bin/bash
#
#SBATCH --partition=ptfgen_uag
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=DAP
#SBATCH --output=../slurm_logs/DAP_consensus_%A.log
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

module load EasyBuild
module load R/4.2.2-foss-2022b

chr=$1 # from "chr1" to "chr22"
analysis=$2 # "DAP_consensus"
dataset=$3 # "blood" or "biopsies"
iteration=$4 # "first" or "second"

ThetaDir=/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/theta/${dataset}/meta_analysis_regulatory_${iteration}
cd ${ThetaDir}
modules=$(find *_${chr}.txt | sed 's/Meta_analysis_for_module_//' | sed 's/_chr.*//' | uniq)
cd "/home/u/u230859/Analysis/eQTL/scripts"

for module in ${modules[@]}
do
	echo ${module}
	Rscript 4_compute_theta.R "DAP_consensus" ${dataset} "CD" ${module} ${chr} "manual" 0 ${iteration}
	Rscript 4_compute_theta.R "DAP_consensus" ${dataset} "UC" ${module} ${chr} "manual" 0 ${iteration}
	Rscript 4_compute_theta.R "DAP_consensus" ${dataset} "IBD" ${module} ${chr} "manual" 0 ${iteration}
done

echo -e "\n"
echo "Computation of theta for colocalization completed"
echo -e "\n"
date
