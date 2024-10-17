#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=concat_vcf
#SBATCH --output=Concat_genotypes_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24000

set -e
set -u
set -o pipefail

date
echo -e "\n"

VCF='2022_11_SYSCID_CROHN_423'
WorkingDir='/massstorage/URT/GEN/UAG/IBD'
InputDir=${WorkingDir}/SYSCID/Helene/Genotypes_QC/5_Impute_genotypes/Output/${VCF}/Dose
OutputDir=${WorkingDir}/SYSCID/Helene/Genotypes_QC/5_Impute_genotypes/Output/${VCF}/Dose_concat
mkdir -p ${OutputDir}

Suffix='.dose.vcf.gz'

# 1) Concatenate the 23 chromosome files
bcftools concat ${InputDir}/chr1${Suffix} ${InputDir}/chr2${Suffix} ${InputDir}/chr3${Suffix} ${InputDir}/chr4${Suffix} ${InputDir}/chr5${Suffix} ${InputDir}/chr6${Suffix} \
${InputDir}/chr7${Suffix} ${InputDir}/chr8${Suffix} ${InputDir}/chr9${Suffix} ${InputDir}/chr10${Suffix} ${InputDir}/chr11${Suffix} ${InputDir}/chr12${Suffix} \
${InputDir}/chr13${Suffix} ${InputDir}/chr14${Suffix} ${InputDir}/chr15${Suffix} ${InputDir}/chr16${Suffix} ${InputDir}/chr17${Suffix} ${InputDir}/chr18${Suffix} \
${InputDir}/chr19${Suffix} ${InputDir}/chr20${Suffix} ${InputDir}/chr21${Suffix} ${InputDir}/chr22${Suffix} ${InputDir}/chrX${Suffix} \
-o ${OutputDir}/CEDAR_concat.vcf --threads ${SLURM_CPUS_PER_TASK} 

bgzip -c ${OutputDir}/CEDAR_concat.vcf > ${OutputDir}/CEDAR_concat.vcf.gz
rm ${OutputDir}/CEDAR_concat.vcf

echo -e "\n"
echo "Concatenation of genotypes completed"
echo -e "\n"
date
