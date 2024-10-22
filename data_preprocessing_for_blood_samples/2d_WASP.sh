#!/bin/bash
#
#SBATCH --partition=all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=WASP
#SBATCH --output=../slurm_logs/3_WASP_%A_%a_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000
#
##SBATCH --array=1-96%8

set -e
set -u
set -o pipefail

date
echo -e "\n"
echo "node used: $SLURM_JOB_NODELIST"
echo "memory per CPU: $SLURM_MEM_PER_CPU"
echo -e "\n"

module load EasyBuild
module load SAMtools/1.17-GCC-12.2.0 #samtools/1.9
module load python/3.7 #module load Python/3.6.6-foss-2018b
#module load jdk/jdk1.8.0_66
module load picard-tools/2.7.1
Picard=$(which picard-tools_2.7.1.jar)

Plate=$1
Tag=$2
Replacement=$3

WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene"
InputDir=${WorkingDir}/RNAseq/data/fastq/${Plate}
AlnDir=${WorkingDir}/RNAseq/aln/${Plate}/STAR_BAM
OutputDir=${WorkingDir}/RNAseq/aln/${Plate}
mkdir -p ${OutputDir}/WASP_BAM
mkdir -p ${OutputDir}/WASP_temp_counts
mkdir -p ${OutputDir}/WASP_counts

ScratchDir="/local/u230859/${SLURM_JOB_ID}/WASP_filter"
mkdir -p ${ScratchDir}

RefGenomeDir="/massstorage/URT/GEN/UAG/IBD/PLATFORMS/GEN/Homo_sapiens/Ensembl/GRCh38/release_105"
GTF=${RefGenomeDir}/Annotation/Genes/Homo_sapiens.GRCh38.105.gtf

VCFDir=${WorkingDir}/Genotypes_QC/5_Impute_genotypes/Output/2022_11_SYSCID_CROHN_423/Dose_concat_filtered
BED_file=${VCFDir}/CEDAR-clean-snps-4_indels.bed

suffix="_R1_001.fastq.gz"
files=( $(find ${InputDir}/*${suffix} -type f -size +10000000c) )
i=$((SLURM_ARRAY_TASK_ID-1))
reads1=${files[$i]}

prefix_temp1=$(basename ${reads1})
prefix_temp2=${prefix_temp1/${suffix}/}
prefix_temp3=${prefix_temp2/${Tag}/${Replacement}}
prefix_temp4=${prefix_temp3/SYSCID_/S-}
prefix_temp5=${prefix_temp4/_Sample_/_}
prefix=${prefix_temp5/S348/S-348}

if [ "$Plate" == 'plate-35' ] & [ "$prefix_temp2" == "S-239_22_AHNG5KDSX5_S464_L002" ]; then
	prefix="S-239bis-22"

elif [ "$Plate" == 'plate-60' ] & [ "$prefix_temp2" == "S-231-8_BHKYVKDSX7_S237_L001" ]; then
	prefix="S-231-8_L1"

elif [ "$Plate" == 'plate-60' ] & [ "$prefix_temp2" == "S-231-8_BHKYVKDSX7_S237_L002" ]; then
	prefix="S-231-8_L2"
fi

echo "Plate: $Plate"
echo "Sample processed: $prefix"
echo -e "\n"

if [ "${prefix}" = "Syscid-218-Maxwell-22-PBMC" ] || [ "${prefix}" = "Syscid-218-Maxwell-23-neutrophils" ] || [ "${prefix}" = "Syscid-218-Maxwell-25-granulocytes" ]
then
    Strandness="reverse"
else
    Strandness="no"
fi

###############################

size=$(du -sb ${reads1} | awk '{print $1}')

if [ "${size}" -ge "4000000000" ]; then
	start_sample=${ScratchDir}/${prefix}_downsampled.bam
	samtools view -h -b -s 1234.1593481 ${AlnDir}/${prefix}_Aligned.sortedByCoord.out.bam > ${start_sample} #| samtools sort -h -b
	samtools index ${start_sample}
	
	echo -e "\n"
	echo "Downsampling completed"
	echo -e "\n"
	date
else
	start_sample=${AlnDir}/${prefix}_Aligned.sortedByCoord.out.bam
fi

###############################

if [ ! -s "${OutputDir}/WASP_counts/${prefix}_ReadsPerGene.out.tab" ]; then
	samtools view ${start_sample} | grep -E "vW:i:2|vW:i:3|vW:i:4|vW:i:5|vW:i:6|vW:i:7" | awk '{print $1}' | sort | uniq > ${ScratchDir}/${prefix}_fail_wasp_filtering.txt || true
	wc -l ${ScratchDir}/${prefix}_fail_wasp_filtering.txt
	samtools view -h ${start_sample} | grep -f ${ScratchDir}/${prefix}_fail_wasp_filtering.txt -F -v -w | samtools view -h -b > ${ScratchDir}/${prefix}_Aligned.sortedByCoord.out.bam
	samtools index ${ScratchDir}/${prefix}_Aligned.sortedByCoord.out.bam
	htseq-count -f bam -r pos -s ${Strandness} -a 10 ${ScratchDir}/${prefix}_Aligned.sortedByCoord.out.bam ${GTF} > ${OutputDir}/WASP_temp_counts/${prefix}_ReadsPerGene.out.tab
	
	###############################
	
	samtools view -M -L ${BED_file} ${ScratchDir}/${prefix}_Aligned.sortedByCoord.out.bam | awk '{print $1}' | sort | uniq > ${ScratchDir}/${prefix}_overlap_indels.txt || true
	wc -l ${ScratchDir}/${prefix}_overlap_indels.txt
	samtools view -h ${ScratchDir}/${prefix}_Aligned.sortedByCoord.out.bam | grep -f ${ScratchDir}/${prefix}_overlap_indels.txt -F -v -w | samtools view -h -b > ${OutputDir}/WASP_BAM/${prefix}_Aligned.sortedByCoord.out.bam
	samtools index ${OutputDir}/WASP_BAM/${prefix}_Aligned.sortedByCoord.out.bam
	htseq-count -f bam -r pos -s ${Strandness} -a 10 ${OutputDir}/WASP_BAM/${prefix}_Aligned.sortedByCoord.out.bam ${GTF} > ${OutputDir}/WASP_counts/${prefix}_ReadsPerGene.out.tab
fi

echo -e "\n"
echo "Removing reads that don't pass WASP filtering and that overlap indels completed"
echo -e "\n"
date
