#!/bin/bash
#
#SBATCH --partition=all_5days,all_10days
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=DAP
#SBATCH --output=../slurm_logs/DAP_%A_%a_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#
#SBATCH --array=1-27
##SBATCH --array=1-401

set -e
set -u
set -o pipefail

date
echo -e "\n"
echo "node used: $SLURM_JOB_NODELIST"

module load EasyBuild
module load R/4.2.2-foss-2022b

chr=$1 # from "chr1" to "chr22"
analysis=$2 # "DAP"
dataset=$3 # "blood" or "biopsies"

B=10000
alpha=0.1
DAP_sym_border=0

echo ${chr}

if [ "$dataset" == "blood" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA"
	
elif [ "$dataset" == "biopsies" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap"
	
fi

ResultsDir=${WorkingDir}/results
cell_types=( ${ResultsDir}/residualised_expression_* )
cell_type_temp1=${cell_types[$((SLURM_ARRAY_TASK_ID-1))]}
cell_type_temp2=$(basename ${cell_type_temp1})
cell_type_temp3=${cell_type_temp2/residualised_expression_/}
cell_type_temp4=${cell_type_temp3/ENSG_/}
cell_type_temp5=${cell_type_temp4/_none_*/}
cell_type=${cell_type_temp5/_full_peer/}
echo ${cell_type}

Rscript 4_compute_theta.R ${analysis} ${dataset} "CD" ${cell_type} ${chr} ${B} ${alpha} "manual" ${DAP_sym_border}
Rscript 4_compute_theta.R ${analysis} ${dataset} "UC" ${cell_type} ${chr} ${B} ${alpha} "manual" ${DAP_sym_border}
Rscript 4_compute_theta.R ${analysis} ${dataset} "IBD" ${cell_type} ${chr} ${B} ${alpha} "manual" ${DAP_sym_border}

echo -e "\n"
echo "Computation of theta for colocalization completed"
echo -e "\n"
date
