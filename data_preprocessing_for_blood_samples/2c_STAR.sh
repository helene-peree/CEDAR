#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=STAR
#SBATCH --output=../slurm_logs/3_STAR_%A_%a_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12000
#
##SBATCH --array=1-96%8

set -e
set -u
set -o pipefail

module load EasyBuild
module load SAMtools/1.17-GCC-12.2.0 #samtools/1.9
module load singularity/3.7.1 #singularity/3.2.1

Plate=$1
Tag=$2
Replacement=$3

seqtk="/home/u/u230859/.local/bin/seqtk/seqtk"
StarSingularityFile='/massstorage/RES/Tools/Singularity/containers/star_v2-7-1a.sif'
WorkingDir="/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene"
InputDir=${WorkingDir}/RNAseq/data/fastq/${Plate}
OutputDir=${WorkingDir}/RNAseq/aln/${Plate}
mkdir -p ${OutputDir}

ScratchDir="/local/u230859/${SLURM_JOB_ID}"
mkdir -p ${ScratchDir}/STAR_BAM
mkdir -p ${ScratchDir}/STAR_counts
mkdir -p ${ScratchDir}/STAR_logs

RefGenomeDir="/massstorage/URT/GEN/UAG/IBD/PLATFORMS/GEN/Homo_sapiens/Ensembl/GRCh38/release_105"
GTF=${RefGenomeDir}/Annotation/Genes/Homo_sapiens.GRCh38.105.gtf
STARIndexDir=${RefGenomeDir}/Sequence/star_v2-7-1a

VCFDir=${WorkingDir}/Genotypes_QC/5_Impute_genotypes/Output/2022_11_SYSCID_CROHN_423/Dose_concat_filtered
VCF_file=CEDAR-clean-snps-4.vcf.gz

export SINGULARITY_BIND=${InputDir},${ScratchDir},${STARIndexDir},${VCFDir}

###############################

date
echo -e "\n"
echo "node: $SLURM_JOB_NODELIST"
echo "job ID: $SLURM_JOB_ID"
echo "number of tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "memory per CPU: $SLURM_MEM_PER_CPU"
echo "array: $SLURM_ARRAY_TASK_ID"
echo -e "\n"

echo "STAR version:"
singularity exec ${StarSingularityFile} STAR --version
module list
echo -e "\n"

###############################

suffix="_R1_001.fastq.gz"
files=( $(find ${InputDir}/*${suffix} -type f -size +10000000c) )
i=$((SLURM_ARRAY_TASK_ID-1))
reads1=${files[$i]}
reads2=${reads1/R1/R2}

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

###############################

size=$(du -sb ${reads1} | awk '{print $1}')

if [ "${size}" -ge "4000000000" ]; then
	start_reads1=${ScratchDir}/${prefix}_downsampled_R1.fastq.gz
	start_reads2=${ScratchDir}/${prefix}_downsampled_R2.fastq.gz
	
	zgrep -H -c @A00 ${reads1}
	zgrep -H -c @A00 ${reads2}
	
	${seqtk} sample -s100 ${reads1} 26184065 | gzip > ${start_reads1}
	${seqtk} sample -s100 ${reads2} 26184065 | gzip > ${start_reads2}
	
	zgrep -H -c @A00 ${start_reads1}
	zgrep -H -c @A00 ${start_reads2}
	
	echo -e "\n"
	echo "Downsampling completed"
	echo -e "\n"
	date
else
	start_reads1=${reads1}
	start_reads2=${reads2}
fi

###############################

if [ ! -f "${VCFDir}/${VCF_file/.vcf.gz/_indels.bed}" ]; then
	bcftools view --types indels --threads ${SLURM_CPUS_PER_TASK} ${VCFDir}/${VCF_file} | awk '/#CHROM/,EOF {print $1"\t"$2-1"\t"$2}' | awk '{sub(/^chr/,"",$0)} 1' | sed 1d > ${VCFDir}/${VCF_file/.vcf.gz/_indels.bed}
	head ${VCFDir}/${VCF_file/.vcf.gz/_indels.bed}
	wc -l ${VCFDir}/${VCF_file/.vcf.gz/_indels.bed}
	
	echo "Extracting 403211 indels completed"
	echo -e "\n"
fi

if [ ! -f "${VCFDir}/${VCF_file/.vcf.gz/_snps.vcf}" ]; then
	cat <(bcftools view -h ${VCFDir}/${VCF_file} | sed '/#CHROM/s/.*/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE/') <(bcftools view -H --types snps --threads ${SLURM_CPUS_PER_TASK} ${VCFDir}/${VCF_file} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t0\/1"}') | awk '{sub(/^chr/,"",$0)} 1' > ${VCFDir}/${VCF_file/.vcf.gz/_snps.vcf}
	head ${VCFDir}/${VCF_file/.vcf.gz/_snps.vcf}
	wc -l ${VCFDir}/${VCF_file/.vcf.gz/_snps.vcf}
	
	echo "Extracting 5896787 SNPs completed"
	echo -e "\n"
fi

###############################

srun singularity exec ${StarSingularityFile} STAR \
--runMode alignReads \
--runThreadN ${SLURM_CPUS_PER_TASK} \
--genomeDir ${STARIndexDir} \
--readFilesIn ${start_reads1} ${start_reads2} \
--readFilesCommand gunzip -c \
--quantMode GeneCounts \
--varVCFfile ${VCFDir}/${VCF_file/.vcf.gz/_snps.vcf} \
--waspOutputMode SAMtag \
--outFileNamePrefix ${ScratchDir}/${prefix}_ \
--outSAMattributes NH HI AS nM vA vG vW \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--limitBAMsortRAM 30000000000

echo -e "\n"
echo "Mapping to genome with STAR completed"
echo -e "\n"
date

###############################

mv ${ScratchDir}/*_Aligned.sortedByCoord.out.bam ${ScratchDir}/STAR_BAM
mv ${ScratchDir}/*_ReadsPerGene.out.tab ${ScratchDir}/STAR_counts
mv ${ScratchDir}/*_Log* ${ScratchDir}/STAR_logs
mv ${ScratchDir}/*_SJ.out.tab* ${ScratchDir}/STAR_logs

if [ "${size}" -ge "4000000000" ]; then
	rm ${start_reads1}
	rm ${start_reads2}
fi

srun rsync -r -p ${ScratchDir}/* ${OutputDir}/.
rm -rf ${ScratchDir}

echo -e "\n"
echo "Moving files to mass storage completed"
echo -e "\n"
date

###############################

samtools index ${OutputDir}/STAR_BAM/${prefix}_Aligned.sortedByCoord.out.bam

echo -e "\n"
echo "Indexing file completed"
echo -e "\n"
date
