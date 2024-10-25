#!/bin/bash
#
#SBATCH --partition=all_5days,all_10days
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=EAP
#SBATCH --output=../slurm_logs/EAP_%A_%a_%j.log
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

chr=$1 # from "chr1" to "chrX"
analysis=$2 # "EAP" or "EAP_manual" or "EAP_consensus"
dataset=$3 # "blood" or "biopsies" or "both"

B=0
alpha=0.1

if [ "$dataset" == "blood" ] || [ "$dataset" == "biopsies" ]; then
	
	if [ "$dataset" == "blood" ]; then
		ResultsDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA/results"
	
	elif [ "$dataset" == "biopsies" ]; then
		ResultsDir="/massstorage/URT/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap/results"
	fi
	
	cell_types=( ${ResultsDir}/residualised_expression_* )
	n=$((${#cell_types[@]}-1))
	i=$((SLURM_ARRAY_TASK_ID-1))
	
	cell_type_temp1=${cell_types[$i]}
	cell_type_temp2=$(basename ${cell_type_temp1})
	cell_type_temp3=${cell_type_temp2/residualised_expression_/}
	cell_type_temp4=${cell_type_temp3/ENSG_/}
	cell_type_temp5=${cell_type_temp4/_none_*/}
	cell_type_1=${cell_type_temp5/_full_peer/}
	cell_types_2=( $(awk -v awkvar=${cell_type_1} '{if ($1 == awkvar) print $2}' balance_number_comparisons_${dataset}.txt) )
	
	for cell_type_2 in ${cell_types_2[@]}
	do
		echo -e "\n"
		echo ${cell_type_1} vs ${cell_type_2}
		Rscript 4_compute_theta.R ${analysis} ${dataset} ${cell_type_1} ${cell_type_2} ${chr} ${B} ${alpha}
	done
	
elif [ "$dataset" == "both" ]; then
	ResultsDir_1="/massstorage/URT/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap/results"
	cell_types_1=( ${ResultsDir_1}/residualised_expression_* )
	cell_type_temp1=${cell_types_1[$((SLURM_ARRAY_TASK_ID-1))]}
	cell_type_temp2=$(basename ${cell_type_temp1})
	cell_type_temp3=${cell_type_temp2/residualised_expression_/}
	cell_type_1=${cell_type_temp3/_full_peer/}
	
	ResultsDir_2="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA/results"
	cell_types_2=( ${ResultsDir_2}/residualised_expression_* )
	
	for cell_type_temp1 in ${cell_types_2[@]}
	do
		cell_type_temp2=$(basename ${cell_type_temp1})
		cell_type_temp3=${cell_type_temp2/residualised_expression_ENSG_/}
		cell_type_2=${cell_type_temp3/_none_*/}
		echo -e "\n"
		echo ${cell_type_1} vs ${cell_type_2}
		Rscript 4_compute_theta.R ${analysis} ${dataset} ${cell_type_1} ${cell_type_2} ${chr} ${B} ${alpha}
	done
fi

echo -e "\n"
echo "Computation of theta between cell types completed"
echo -e "\n"
date
