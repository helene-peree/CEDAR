#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=ruler
#SBATCH --output=../slurm_logs/0_ruler_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

set -e
set -u
set -o pipefail

date
echo -e "\n"

Plate=$1
Nb_samples_at_once=$2
InputDir=/massstorage/URT/GEN/UAG/IBD/SYSCID/Helene/RNAseq/data/fastq/${Plate}

##################################################

Nb_samples_start=1
Nb_samples_end=$(find ${InputDir}/*_R1_001.fastq.gz -type f -size +10000000c | wc -l)
Nb_samples_total=$(find ${InputDir}/*_R1_001.fastq.gz -type f | wc -l)

echo "Plate: $Plate"
echo "Number of samples processed: $Nb_samples_end"
echo "Number of samples sequenced: $Nb_samples_total"
echo -e "\n"

if [ "$Plate" == 'plate-A' ] || [ "$Plate" == 'plate-D' ] || [ "$Plate" == 'plate-E' ] || [ "$Plate" == 'plate-J' ]; then
	Tag="_NGS*_L00"
	Replacement="_L"
	
elif [ "$Plate" == 'plate-B' ] || [ "$Plate" == "plate-C" ] || [ "$Plate" == 'plate-F' ] || [ "$Plate" == 'plate-G' ] || [ "$Plate" == 'plate-H' ] || [ "$Plate" == 'plate-I' ]; then
	Tag="_NGS*"
	Replacement=""
	
elif [ "$Plate" == 'plate-0' ]; then
	Tag="_NGS20*"
	Replacement=""
	
elif [ "$Plate" == 'plate-1' ] || [ "$Plate" == 'plate-2' ]; then
	Tag="_AHJVHJDSX2*"
	Replacement=""
	
elif [ "$Plate" == 'plate-3' ] || [ "$Plate" == 'plate-4' ] || [ "$Plate" == 'plate-5' ] || [ "$Plate" == 'plate-6' ]; then
	Tag="_AHVVFLDSX2*"
	Replacement=""
	
elif [ "$Plate" == 'plate-7' ] || [ "$Plate" == 'plate-8' ]; then
	Tag="_AH25VNDSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-9' ] || [ "$Plate" == 'plate-10' ] || [ "$Plate" == 'plate-11' ] || [ "$Plate" == 'plate-12' ]; then
	Tag="_AH5TGJDSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-13' ]; then
	Tag="_AH7WM3DSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-14' ] || [ "$Plate" == 'plate-15' ] || [ "$Plate" == 'plate-16' ]; then
	Tag="_AH7TV2DSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-17' ] || [ "$Plate" == 'plate-18' ]; then
	Tag="_AHJH7CDSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-19' ] || [ "$Plate" == 'plate-20' ]; then
	Tag="-PL*"
	Replacement=""
	
elif [ "$Plate" == 'plate-21' ] || [ "$Plate" == 'plate-22' ] || [ "$Plate" == 'plate-23' ] || [ "$Plate" == 'plate-25' ]; then
	Tag="_BHJGTNDSX3*"
	Replacement=""
	
elif [ "$Plate" == 'plate-27' ] || [ "$Plate" == 'plate-29' ] || [ "$Plate" == 'plate-31' ]; then
	Tag="_BH3MM2DSX5*"
	Replacement=""
	
elif [ "$Plate" == 'plate-33' ] || [ "$Plate" == 'plate-34' ] || [ "$Plate" == 'plate-37' ] || [ "$Plate" == 'plate-38' ]; then
	Tag="_AHCGWLDSX5*"
	Replacement=""
	
elif [ "$Plate" == 'plate-35' ] || [ "$Plate" == 'plate-39' ] || [ "$Plate" == 'plate-40' ]; then
	Tag="_AHNG5KDSX5*"
	Replacement=""
	
elif [ "$Plate" == 'plate-41' ] || [ "$Plate" == 'plate-43' ] || [ "$Plate" == 'plate-45' ]; then
	Tag="_AHVNKTDSX5*"
	Replacement=""
	
elif [ "$Plate" == 'plate-47' ] || [ "$Plate" == 'plate-49' ]; then
	Tag="_AH27KHDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-51' ] || [ "$Plate" == 'plate-52' ] || [ "$Plate" == 'plate-53' ] || [ "$Plate" == 'plate-54' ]; then
	Tag="_AH2NWMDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-55' ] || [ "$Plate" == 'plate-56' ] || [ "$Plate" == 'plate-57' ] || [ "$Plate" == 'plate-58' ]; then
	Tag="_BHHC3MDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-59' ] || [ "$Plate" == 'plate-60' ] || [ "$Plate" == 'plate-61' ] || [ "$Plate" == 'plate-62' ]; then
	Tag="_BHKYVKDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-63' ] || [ "$Plate" == 'plate-64' ] || [ "$Plate" == 'plate-67' ] || [ "$Plate" == 'plate-68' ]; then
	Tag="_AHJ52YDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-65' ] || [ "$Plate" == 'plate-66' ] || [ "$Plate" == 'plate-75' ] || [ "$Plate" == 'plate-76' ]; then
	Tag="_AH5G5HDSXC*"
	Replacement=""
	
elif [ "$Plate" == 'plate-71' ] || [ "$Plate" == 'plate-72' ]; then
	Tag="_AHCY7CDSXC*"
	Replacement=""
	
elif [ "$Plate" == 'plate-73' ] || [ "$Plate" == 'plate-74' ] || [ "$Plate" == 'plate-77' ] || [ "$Plate" == 'plate-78' ]; then
	Tag="_AHHWCNDSX7*"
	Replacement=""
	
elif [ "$Plate" == 'plate-79' ] || [ "$Plate" == 'plate-80' ]; then
	Tag="_AHF77LDSXC*"
	Replacement=""
	
elif [ "$Plate" == 'plate-83' ] || [ "$Plate" == 'plate-85' ]; then
	Tag="_AHHW2JDSXC*"
	Replacement=""
	
elif [ "$Plate" == 'plate-69' ] || [ "$Plate" == 'plate-81' ] || [ "$Plate" == 'plate-87' ]; then
	Tag="_AHF7FWDSXC*"
	Replacement=""
fi

##################################################

echo FASTQC
FASTQC=$(sbatch --array="$Nb_samples_start"-"$Nb_samples_end"%"$Nb_samples_at_once" 2a_fastQC.sh "$Plate" "$Tag" "$Replacement" | awk '{print $4}')
echo $FASTQC

echo MULTIQC_FASTQC
MULTIQC_FASTQC=$(sbatch --dependency=afterok:"$FASTQC" 2b_multiQC_per_plate.sh "$Plate" | awk '{print $4}')
echo $MULTIQC_FASTQC

echo STAR
STAR=$(sbatch --dependency=afterany:"$FASTQC" --array="$Nb_samples_start"-"$Nb_samples_end"%"$Nb_samples_at_once" 2c_STAR.sh "$Plate" "$Tag" "$Replacement" | awk '{print $4}')
echo $STAR

echo PICARD_STAR
PICARD_STAR=$(sbatch --dependency=afterany:"$STAR" --array="$Nb_samples_start"-"$Nb_samples_end"%"$Nb_samples_at_once" 2e_picard.sh "$Plate" "$Tag" "$Replacement" "STAR" | awk '{print $4}')
echo $PICARD_STAR

echo MULTIQC_STAR
MULTIQC_STAR=$(sbatch --dependency=afterok:"$PICARD_STAR" 2f_multiQC_per_plate.sh "$Plate" "STAR" | awk '{print $4}')
echo $MULTIQC_STAR

echo WASP
WASP=$(sbatch --dependency=afterany:"$PICARD_STAR" --array="$Nb_samples_start"-"$Nb_samples_end"%"$Nb_samples_at_once" 2d_WASP.sh "$Plate" "$Tag" "$Replacement" | awk '{print $4}')
echo $WASP

echo PICARD_WASP
PICARD_WASP=$(sbatch --dependency=afterany:"$WASP" --array="$Nb_samples_start"-"$Nb_samples_end"%"$Nb_samples_at_once" 2e_picard.sh "$Plate" "$Tag" "$Replacement" "WASP" | awk '{print $4}')
echo $PICARD_WASP

echo MULTIQC_WASP
MULTIQC_WASP=$(sbatch --dependency=afterok:"$PICARD_WASP" 2f_multiQC_per_plate.sh "$Plate" "WASP" | awk '{print $4}')
echo $MULTIQC_WASP

echo -e "\n"
echo "Scripts launched"
echo -e "\n"
date
