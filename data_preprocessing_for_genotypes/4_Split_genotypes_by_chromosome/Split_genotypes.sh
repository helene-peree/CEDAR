#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=split_vcf
#SBATCH --output=Split_genotypes_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000

set -e
set -u
set -o pipefail

date
echo -e "\n"

module load tabix/0.2.6

VCF='2022_11_SYSCID_CROHN_423'
WorkingDir='/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/Genotypes_QC'
OutputDir=${WorkingDir}/4_Split_genotypes_by_chromosome/Output/${VCF}
mkdir -p ${OutputDir}

InputDir=${WorkingDir}/3_Liftover_genotypes/Output/${VCF}
InputFile=${InputDir}/CEDAR-clean-snps-2_mixed_build.vcf.gz
TempFile=${OutputDir}/CEDAR-clean-snps-2_mixed_build.vcf.bgz

gunzip -c ${InputFile} | bgzip -c > ${TempFile}
tabix -p vcf ${TempFile}

for chromosome in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
	tabix -h ${TempFile} ${chromosome} | bgzip -c > ${TempFile/.vcf.bgz/_${chromosome}.vcf.gz}
done

echo -e "\n"
echo "Number of variants in chrY: "
gunzip -c ${TempFile/.vcf.bgz/_chrY.vcf.gz} | awk '{if($0 !~ /^#/) print $0}' | wc -l
echo -e "\n"
echo "Splitting of genotypes completed"
echo -e "\n"
date
