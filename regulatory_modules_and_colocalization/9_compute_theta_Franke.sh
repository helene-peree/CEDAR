#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days,all_10days
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=Franke
#SBATCH --output=../slurm_logs/Franke_%A_%a_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#
#SBATCH --array=1-22%22

set -e
set -u
set -o pipefail

date
echo -e "\n"
echo "node used: $SLURM_JOB_NODELIST"

module load EasyBuild
module load R/4.2.2-foss-2022b

chrs=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" )
i=$((SLURM_ARRAY_TASK_ID-1))
chr=${chrs[$i]}

DAP_borders="manual"
DAP_sym_border=0

Rscript 9_compute_theta_Franke.R "DAP" "CD" ${chr} ${DAP_borders} ${DAP_sym_border}
Rscript 9_compute_theta_Franke.R "DAP" "IBD" ${chr} ${DAP_borders} ${DAP_sym_border}
Rscript 9_compute_theta_Franke.R "DAP" "UC" ${chr} ${DAP_borders} ${DAP_sym_border}

echo -e "\n"
echo "Computation of theta completed"
echo -e "\n"
date
