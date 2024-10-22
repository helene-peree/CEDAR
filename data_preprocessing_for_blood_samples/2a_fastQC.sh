#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=fastQC
#SBATCH --output=../slurm_logs/1_fastQC_%A_%a_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#
##SBATCH --array=1-96%8

set -e
set -u
set -o pipefail

date
echo -e "\n"

module load fastqc/0.12.1

Plate=$1
Tag=$2
Replacement=$3

WorkingDir='/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene'
InputDir=${WorkingDir}/RNAseq/data/fastq/${Plate}
OutputDir=${WorkingDir}/RNAseq/data/fastQC/${Plate}
mkdir -p ${OutputDir}

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

cd ${OutputDir}

ends=( "R1" "R2" )
for end in "${ends[@]}"
do
	old_name=${prefix_temp2}_${end}_001
	new_name=${prefix}_${end}
	
	if [ ! -s "${new_name}_fastqc.zip" ]; then
		fastqc --extract --nogroup --outdir=${OutputDir} ${InputDir}/${old_name}.fastq.gz
		echo -e "\n"
		
		rm ${old_name}_fastqc.zip
		mv ${old_name}_fastqc ${new_name}_fastqc
		mv ${old_name}_fastqc.html ${new_name}_fastqc.html
		
		sed -i "s|${old_name}.fastq.gz|${new_name}|g" ${new_name}_fastqc/fastqc.fo
		sed -i "s|${old_name}.fastq.gz|${new_name}|g" ${new_name}_fastqc/fastqc_data.txt
		sed -i "s|${old_name}.fastq.gz|${new_name}|g" ${new_name}_fastqc/fastqc_report.html
		sed -i "s|${old_name}.fastq.gz|${new_name}|g" ${new_name}_fastqc/summary.txt
		sed -i "s|${old_name}.fastq.gz|${new_name}|g" ${new_name}_fastqc.html
		
		zip -r ${new_name}_fastqc.zip ${new_name}_fastqc
		rm -r ${new_name}_fastqc
		echo -e "\n"
	fi
done

echo -e "\n"
echo "QC analysis with fastQC completed"
echo -e "\n"
date
