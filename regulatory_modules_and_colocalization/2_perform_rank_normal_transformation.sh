#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days,all_10days
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=normal
#SBATCH --output=../slurm_logs/normal_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=56000

set -e
set -u
set -o pipefail

date
echo -e "\n"
echo "node used: $SLURM_JOB_NODELIST"

module load qtltools

dataset=$1 # "blood" or "biopsies"

if [ "$dataset" == "blood" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA"
	suffix="_none_*"
	replacement="_none"
	
elif [ "$dataset" == "biopsies" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap"
	suffix="_full_peer"
	replacement=""
fi

ResultsDir=${WorkingDir}/results
cell_types=( ${ResultsDir}/residualised_expression_* )

for cell_type_temp1 in ${cell_types[@]}
do
	cell_type_temp2=$(basename ${cell_type_temp1})
	cell_type_temp3=${cell_type_temp2/residualised_expression_/}
	cell_type=${cell_type_temp3/${suffix}/${replacement}}
	echo ${cell_type}
	
	ExpressionDir=${ResultsDir}/${cell_type_temp2}
	
	if [ ! -f ${ExpressionDir}/expression_rank_normal_transformed_${cell_type}.bed ]; then
		QTLtools correct --bed ${ExpressionDir}/expression_residualised_${cell_type} --normal --out ${ExpressionDir}/expression_rank_normal_transformed_${cell_type}.bed
	fi
done

echo -e "\n"
echo "Performing rank normal transformation for expression data completed"
echo -e "\n"
date
