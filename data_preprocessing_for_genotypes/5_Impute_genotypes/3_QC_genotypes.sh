#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=QC_vcf
#SBATCH --output=QC_genotypes_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32000

set -u
set -o pipefail

date
echo -e "\n"

module load R
module load plink/1.9
module load qtltools/1.3.1

VCFDir='2022_11_SYSCID_CROHN_423'
Threshold_r2=0.7
Threshold_maf=0.05

MetaDir='/home/u/u230859/Analysis/Clinical_Fortessa_data'
ScriptDir='/home/u/u230859/Analysis/Genotypes_QC/5_Impute_genotypes'
WorkingDir=/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene
InputDir=${WorkingDir}/Genotypes_QC/5_Impute_genotypes/Output/${VCFDir}/Dose_concat
OutputDir=${WorkingDir}/Genotypes_QC/5_Impute_genotypes/Output/${VCFDir}/Dose_concat_filtered
mkdir -p ${OutputDir}

Sex=${MetaDir}/Sex_information.txt
Affection=${MetaDir}/Affection_status.txt
VCF=${InputDir}/CEDAR_concat.vcf.gz

Prefix='CEDAR'
FileCleaned1=${OutputDir}/${Prefix}-clean-snps-1
FileCleaned2=${OutputDir}/${Prefix}-clean-snps-2
FileCleaned3=${OutputDir}/${Prefix}-clean-snps-3
FileCleaned4=${OutputDir}/${Prefix}-clean-snps-4

# 1) Remove CROHN and COV samples
# 2) Remove FORMAT fields other than GT and DS
if [ "$VCFDir" == "2022_11_SYSCID_60" ]; then
	bcftools view -S ${MetaDir}/Biopsies_order.txt -Oz --threads ${SLURM_CPUS_PER_TASK} ${VCF} > ${FileCleaned1}.vcf.gz
	bcftools annotate -x ^FORMAT/GT,FORMAT/DS -Oz --threads ${SLURM_CPUS_PER_TASK} ${FileCleaned1}.vcf.gz > ${FileCleaned2}.vcf.gz
	
elif [ "$VCFDir" == "2022_11_SYSCID_CROHN_423" ]; then
	zgrep -F -w -m 1 "#CHROM" ${VCF} | tr "\t" "\n" | grep -E "CROHN|COV" > ${OutputDir}/crohn_cov.txt
	bcftools view -S ^${OutputDir}/crohn_cov.txt -Oz --threads ${SLURM_CPUS_PER_TASK} ${VCF} > ${FileCleaned1}.vcf.gz
	bcftools annotate -x ^FORMAT/GT,FORMAT/DS -Oz --threads ${SLURM_CPUS_PER_TASK} ${FileCleaned1}.vcf.gz > ${FileCleaned2}.vcf.gz
fi

# 3) Remove markers with imputation R2 ≤ 0.7
bcftools view -i "R2>${Threshold_r2}" -Oz --threads ${SLURM_CPUS_PER_TASK} ${FileCleaned2}.vcf.gz > ${FileCleaned3}.vcf.gz

# 4) Create BED, BIM and FAM files
plink --vcf ${FileCleaned3}.vcf.gz --update-sex ${Sex} --pheno ${Affection} --keep-allele-order --make-bed --out ${FileCleaned3}

# 5) Remove markers with MAF < 0.05
plink --bfile ${FileCleaned3} --freq --maf ${Threshold_maf} --keep-allele-order --make-bed --out ${FileCleaned4}
Rscript ${ScriptDir}/3_MAF.R ${FileCleaned4} ${Threshold_maf}
zgrep -f ${OutputDir}/fail-qc-maf.txt -F -v -w ${FileCleaned3}.vcf.gz > ${FileCleaned4}.vcf

# 6) Prepare files for nextflow eQTL pipeline
bgzip -c ${FileCleaned4}.vcf > ${FileCleaned4}.vcf.gz
tabix -p vcf ${FileCleaned4}.vcf.gz
QTLtools pca --vcf ${FileCleaned4}.vcf.gz --scale --center --out ${OutputDir}/genotypes_PCA

echo -e "\n"
echo "Quality check of genotypes completed"
echo -e "\n"
date
