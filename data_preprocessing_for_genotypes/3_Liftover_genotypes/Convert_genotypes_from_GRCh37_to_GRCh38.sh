#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=convert_vcf
#SBATCH --output=Convert_genotypes_from_GRCh37_to_GRCh38_%A.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16000

set -e
set -u
set -o pipefail

date
echo -e "\n"

module load jdk/jdk1.8.0_66
module load picard-tools/2.7.1
Picard=`which picard-tools_2.7.1.jar`

VCF='2022_11_SYSCID_CROHN_423'
WorkingDir='/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/Genotypes_QC'
OutputDir=${WorkingDir}/3_Liftover_genotypes/Output/${VCF}
mkdir -p ${OutputDir}

Number="2"
CHAINF=/home/u/u230859/Analysis/Genotypes_QC/3_Liftover_genotypes/GRCh37_to_GRCh38.chain.gz
INPUT=${WorkingDir}/2_QC_genotypes/Output/${VCF}/CEDAR-clean-snps-${Number}.vcf.gz
REF='/massstorage/URT/GEN/UAG/IBD/PLATFORMS/GEN/Homo_sapiens/Ensembl/GRCh38/release_105/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa' 
OUTPUT1=${OutputDir}/CEDAR-clean-snps-${Number}_GRCh38.vcf
OUTPUT2=${OutputDir}/CEDAR-clean-snps-${Number}_mixed_build.vcf
RVCF=${OutputDir}/CEDAR-clean-snps-${Number}_rejected.vcf

java -jar ${Picard} LiftoverVcf \
     I=${INPUT} \
	 CHAIN=${CHAINF} \
     O=${OUTPUT1}.gz \
     REJECT=${RVCF}.gz \
     R=${REF}

gunzip -c ${OUTPUT1}.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > ${OUTPUT2}
bgzip -c ${OUTPUT2} > ${OUTPUT2}.gz

echo -e "\n"
echo "Number of variants not lifted over: "
gunzip -c ${RVCF}.gz | awk '{if($0 !~ /^#/) print $0}' | wc -l
echo -e "\n"
echo "Lifting over of genotypes completed"
echo -e "\n"
date
