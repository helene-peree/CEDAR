#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=QC_vcf
#SBATCH --output=QC_genotypes_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000

set -u
set -o pipefail

date
echo -e "\n"

module load plink/1.9
module load R

Batch_VCF='2022_11_SYSCID_CROHN_423'

Threshold_missing_inds=0.03
Threshold_missing_snps=0.05
Threshold_hwe=0.0003

HomeDir='/home/u/u230859/Analysis'
MassDir='/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene'
MetaDir=${HomeDir}/Clinical_Fortessa_data
ScriptDir=${HomeDir}/Genotypes_QC/2_QC_genotypes
DataDir=${MassDir}/Genotypes_QC/1_Merge_genotypes/Output/${Batch_VCF}
HapMapDir=${MassDir}/Genotypes_QC/2_QC_genotypes/HapMap_article
OutputDir=${MassDir}/Genotypes_QC/2_QC_genotypes/Output/${Batch_VCF}

VCF=${DataDir}/Genotypes_merged_${Batch_VCF}
Sex=${MetaDir}/Sex_information.txt
Affection=${MetaDir}/Affection_status.txt

Prefix='CEDAR'
FileRaw=${OutputDir}/${Prefix}-raw
FileNoMissing=${OutputDir}/${Prefix}-no-missing
FileFiltered=${OutputDir}/${Prefix}-filtered
FileMerged=${OutputDir}/${Prefix}-hapmap3
FileCleanedInds=${OutputDir}/${Prefix}-clean-inds
FileCleaned1=${OutputDir}/${Prefix}-clean-snps-1
FileCleaned2=${OutputDir}/${Prefix}-clean-snps-2

cp ${MetaDir}/fail-qc-inclusion.txt ${OutputDir}

# 3) Create BED, BIM and FAM files
plink --vcf ${VCF}_with_header.vcf.gz --keep-allele-order --make-bed --update-sex ${Sex} --pheno ${Affection} --out ${FileRaw}

# 7) Calculate the number and the proportion of missing SNPs per individual (--missing)
# 8) Calculate the observed number of homozygous genotypes and the number of nonmissing genotypes per individual (--het)
# 9) Calculate the heterozygosity rate and create a plot
plink --bfile ${FileRaw} --missing --het --out ${FileRaw}
Rscript ${ScriptDir}/1_imiss-het.R ${FileRaw} ${Threshold_missing_inds}

if [[ ! -e ${OutputDir}/fail-qc-imiss.txt ]]; then
    touch ${OutputDir}/fail-qc-imiss.txt
fi

# 10) Remove individuals with elevated missing (or extreme heterozygosity)
plink --bfile ${FileRaw} --remove ${OutputDir}/fail-qc-imiss.txt --keep-allele-order --make-bed --out ${FileNoMissing}

# 4) Calculate the mean homozygosity rate across X-chromosome markers per individual (--check-sex)
# 5) Identify individuals with discordant sex information
plink --bfile ${FileNoMissing} --check-sex --missing --out ${FileNoMissing}
Rscript ${ScriptDir}/2_sex.R ${FileNoMissing}

# 11) Identify the pairs of SNPs (within a given number of base pairs) that have an r2 value > 0.2
# 12) Calculate pairwise IBD on a subset of markers for all pairs of individuals
# 13) Identify one individual from all pairs of duplicated or related individuals
plink --bfile ${FileNoMissing} --indep-pairwise 50 5 0.2 --out ${FileNoMissing}
plink --bfile ${FileNoMissing} --extract ${FileNoMissing}.prune.in --genome --out ${FileNoMissing}
Rscript ${ScriptDir}/3_IBD.R ${FileNoMissing}

SNPsEasyToAlign=${HapMapDir}/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt
HapMap_CEU_CHB_JPT_YRI=${HapMapDir}/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps

# 14) Keep markers that are not difficult to align
# 15) Merge study genotypes with HapMap data for a subset of markers
plink --bfile ${FileNoMissing} --allow-no-sex --extract ${SNPsEasyToAlign} --keep-allele-order --make-bed --out ${FileFiltered}
plink --bfile ${FileFiltered} --bmerge ${HapMap_CEU_CHB_JPT_YRI} --allow-no-sex --extract ${FileNoMissing}.prune.in --keep-allele-order --make-bed --out ${FileMerged}

if [ -f ${FileMerged}-merge.missnp ]
    then
	# 14) Again with --flip
	# 15) Again
	echo "Flipping strand errors"
	plink --bfile ${FileNoMissing} --allow-no-sex --extract ${SNPsEasyToAlign} --flip ${FileMerged}-merge.missnp --keep-allele-order --make-bed --out ${FileFiltered}
	plink --bfile ${FileFiltered} --bmerge ${HapMap_CEU_CHB_JPT_YRI} --allow-no-sex --extract ${FileNoMissing}.prune.in --keep-allele-order --make-bed --out ${FileMerged}
fi

# 17) Conduct a PCA on the merged data
plink --bfile ${FileMerged} --pca 3 tabs --out ${FileMerged}

# 18) Create a plot of the first 2 PCs
# 19) Exclude individuals who do not match the European population
Rscript ${ScriptDir}/4_PCA.R ${FileMerged}

# 20) List individuals who fail the different QC steps
cat ${OutputDir}/fail-qc-* | sort -k1 | uniq > ${OutputDir}/fail-qc-inds.txt

# 21) Remove all individuals failing QC
plink --bfile ${FileNoMissing} --remove ${OutputDir}/fail-qc-inds.txt --keep-allele-order --make-bed --out ${FileCleanedInds}

# 22) Calculate the missing genotype rate and the HWE p-values
# 26) Remove markers with missing data > 0.05 and HWE P-value < 0.0003
plink --bfile ${FileCleanedInds} --missing --geno ${Threshold_missing_snps} --keep-allele-order --make-bed --out ${FileCleaned1}
plink --bfile ${FileCleaned1} --hardy --hwe ${Threshold_hwe} --keep-allele-order --make-bed --out ${FileCleaned2}

# 23) Plot the missing genotype rate and the HWE
Rscript ${ScriptDir}/6_lmiss.R ${FileCleaned1} ${Threshold_missing_snps}
Rscript ${ScriptDir}/7_HWE.R ${FileCleaned2} ${Threshold_hwe}

# Export VCF file
plink --bfile ${FileCleaned2} --recode vcf-iid --keep-allele-order --output-chr M --out ${FileCleaned2}
gzip ${FileCleaned2}.vcf

echo -e "\n"
echo "Quality check of genotypes completed"
echo -e "\n"
date
