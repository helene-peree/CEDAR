#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=picard
#SBATCH --output=../slurm_logs/4_picard_%A_%a_%j.log
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

#module load jdk/jdk1.8.0_66
module load picard-tools/2.7.1
Picard=$(which picard-tools_2.7.1.jar)

Plate=$1
Tag=$2
Replacement=$3
Software=$4 #STAR or WASP

WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene"
InputDir=${WorkingDir}/RNAseq/data/fastq/${Plate}
AlnDir=${WorkingDir}/RNAseq/aln/${Plate}/${Software}_BAM
OutputDir=${WorkingDir}/RNAseq/aln/${Plate}/${Software}_QC
mkdir -p ${OutputDir}

RefGenomeDir="/massstorage/URT/GEN/UAG/IBD/PLATFORMS/GEN"
FastaRefGen=${RefGenomeDir}/Homo_sapiens/Ensembl/GRCh38/release_105/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa
RefFlat=${RefGenomeDir}/Homo_sapiens/Ensembl/GRCh38/release_105/Annotation/Genes/Homo_sapiens.GRCh38.105_RefFlat.txt

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
    Strandness="SECOND_READ_TRANSCRIPTION_STRAND"
else
    Strandness="NONE"
fi

###############################

srun java -jar ${Picard} CollectAlignmentSummaryMetrics \
           R=${FastaRefGen} \
           I=${AlnDir}/${prefix}_Aligned.sortedByCoord.out.bam \
           O=${OutputDir}/${prefix}_alignment_metrics.out

echo -e "\n"
echo "Collection of alignment metrics with picard completed"
echo -e "\n"
date

###############################

if [ ! -f "${OutputDir}/${prefix}_RNAseq_metrics.out" ]; then
	srun java -jar ${Picard} CollectRnaSeqMetrics \
		   I=${AlnDir}/${prefix}_Aligned.sortedByCoord.out.bam \
		   O=${OutputDir}/${prefix}_RNAseq_metrics.out \
		   REF_FLAT=${RefFlat} \
		   STRAND=${Strandness}
fi

echo -e "\n"
echo "Collection of RNAseq metrics with picard completed"
echo -e "\n"
date
