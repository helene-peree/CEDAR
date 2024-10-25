#!/bin/bash
#
#SBATCH --partition=all_5days,all_10days
#SBATCH --constraint="intel"
#SBATCH --exclude=ptfgen016
#
#SBATCH --job-name=split
#SBATCH --output=../slurm_logs/split_%A.log
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

module load EasyBuild
module load R/4.2.2-foss-2022b

dataset=$1 # "blood" or "biopsies"

if [ "$dataset" == "blood" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA"
	RawFile="CEDAR-clean-snps-4"
	
elif [ "$dataset" == "biopsies" ]; then
	WorkingDir="/massstorage/URT/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap"
	RawFile="CEDAR_sc"
	
fi

RawDir=${WorkingDir}/raw_data

if [ ! -f ${RawDir}/${RawFile}-chr1.rds ]; then
	Rscript 3_split_genotypes_by_chromosome.R ${RawDir} ${RawFile}
fi

echo -e "\n"
echo "Splitting genotypes by chromosome completed"
echo -e "\n"
date
