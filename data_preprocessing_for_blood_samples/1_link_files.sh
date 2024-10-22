#!/bin/bash
#
#SBATCH --partition=all_5hrs,all_24hrs,all_5days
#SBATCH --constraint="intel"
#
#SBATCH --job-name=linking
#SBATCH --output=../slurm_logs/0_link_files_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

set -e
set -u
set -o pipefail

date
echo -e "\n"

WorkingDir='/massstorage/URT/GEN/UAG/IBD'
InputDir_Geneva=${WorkingDir}/PLATFORMS/GEN/RawData_Geneva/downloads/RNA_Dermitzakis_SYSCID_SMARTSeq_Mar21
InputDir_Liege=${WorkingDir}/PLATFORMS/GEN/RawData
OutputDir=${WorkingDir}/SYSCID/Helene/RNAseq/data/fastq
mkdir -p ${OutputDir}

ln -s ${InputDir_Geneva}/211011_A00485_0186_BHLFYHDRXY ${OutputDir}/plate-A
ln -s ${InputDir_Geneva}/211208_A00485_0212_AHLWG2DRXY ${OutputDir}/plate-B
ln -s ${InputDir_Geneva}/211213_A00485_0216_AHLW5NDRXY ${OutputDir}/plate-C
ln -s ${InputDir_Geneva}/220204_A00485_0234_BH5L3NDMXY ${OutputDir}/plate-D
ln -s ${InputDir_Geneva}/220208_A00485_0236_BH5K2TDMXY ${OutputDir}/plate-E
ln -s ${InputDir_Geneva}/220211_A00485_0239_BH7LKJDMXY ${OutputDir}/plate-F
ln -s ${InputDir_Geneva}/220224_A00485_0244_AH7JMVDMXY ${OutputDir}/plate-G
ln -s ${InputDir_Geneva}/220301_A00485_0247_AH7NN5DMXY ${OutputDir}/plate-H
ln -s ${InputDir_Geneva}/220324_A00485_0259_AHTLM3DRXY ${OutputDir}/plate-I
ln -s ${InputDir_Geneva}/220413_A00485_0268_BHH2GVDSX3 ${OutputDir}/plate-J

ln -s ${InputDir_Liege}/200526_A00801_0062_BHFFT3DSXY_20200525-HPeree/fastq ${OutputDir}/plate-0
ln -s ${InputDir_Liege}/210825_A00801_0141_AHJVHJDSX2_20210825-GEN579-Hperee ${OutputDir}/plate-1
ln -s ${InputDir_Liege}/210825_A00801_0141_AHJVHJDSX2_20210825-GEN579-Hperee-Plate2 ${OutputDir}/plate-2
ln -s ${InputDir_Liege}/211119_A00801_0167_AHVVFLDSX2_20211119-GEN0646-Hperee-P3/fastq_bulk ${OutputDir}/plate-3
ln -s ${InputDir_Liege}/211119_A00801_0167_AHVVFLDSX2_20211119-GEN0646-Hperee-P4/fastq_bulk ${OutputDir}/plate-4
ln -s ${InputDir_Liege}/211119_A00801_0167_AHVVFLDSX2_20211119-GEN0646-Hperee-P5/fastq_bulk ${OutputDir}/plate-5
ln -s ${InputDir_Liege}/211119_A00801_0167_AHVVFLDSX2_20211119-GEN0646-Hperee-P6/fastq_bulk ${OutputDir}/plate-6
ln -s ${InputDir_Liege}/211221_A00801_0177_AH25VNDSX3_20211215-GEN646-Hperee-Plate7/fastq_bulk ${OutputDir}/plate-7
ln -s ${InputDir_Liege}/211221_A00801_0177_AH25VNDSX3_20211215-GEN646-Hperee-Plate8/fastq_bulk ${OutputDir}/plate-8
ln -s ${InputDir_Liege}/220113_A00801_0182_AH5TGJDSX3_20220113-GEN646-Hperee-Plate9/fastq_bulk ${OutputDir}/plate-9
ln -s ${InputDir_Liege}/220113_A00801_0182_AH5TGJDSX3_20220113-GEN646-Hperee-Plate10/fastq_bulk ${OutputDir}/plate-10
ln -s ${InputDir_Liege}/220113_A00801_0182_AH5TGJDSX3_20220113-GEN646-Hperee-Plate11/fastq_bulk ${OutputDir}/plate-11
ln -s ${InputDir_Liege}/220113_A00801_0182_AH5TGJDSX3_20220113-GEN646-Hperee-Plate12/fastq_bulk ${OutputDir}/plate-12
ln -s ${InputDir_Liege}/220224_A00801_0200_AH7WM3DSX3_20220217-GEN0809-Hperee-P13/fastq_bulk_All ${OutputDir}/plate-13
ln -s ${InputDir_Liege}/220401_A00801_0210_AH7TV2DSX3_20220331-GEN0880-Hperee-P14/fastq_bulk_All ${OutputDir}/plate-14
ln -s ${InputDir_Liege}/220401_A00801_0210_AH7TV2DSX3_20220331-GEN0880-Hperee-P15/fastq_bulk_All ${OutputDir}/plate-15
ln -s ${InputDir_Liege}/220401_A00801_0210_AH7TV2DSX3_20220331-GEN0880-Hperee-P16/fastq_bulk_All ${OutputDir}/plate-16
ln -s ${InputDir_Liege}/220512_A00801_0224_AHJH7CDSX3_20220509-GEN880-Hperee-Pl17/fastq_bulk_All ${OutputDir}/plate-17
ln -s ${InputDir_Liege}/220512_A00801_0224_AHJH7CDSX3_20220509-GEN880-Hperee-Pl18/fastq_bulk_All ${OutputDir}/plate-18
ln -s ${InputDir_Liege}/220617_A00801_0233_AHJ5KFDSX3_20220615-HPeree-PL19/fastq_bulk_All ${OutputDir}/plate-19
ln -s ${InputDir_Liege}/220617_A00801_0233_AHJ5KFDSX3_20220615-Hperee-PL20/fastq_bulk_All ${OutputDir}/plate-20
ln -s ${InputDir_Liege}/220816_A00801_0259_BHJGTNDSX3_20220811-GEN911-Hperee-Plate21/fastq_bulk_All ${OutputDir}/plate-21
ln -s ${InputDir_Liege}/220816_A00801_0259_BHJGTNDSX3_20220811-GEN911-Hperee-Plate22/fastq_bulk_All ${OutputDir}/plate-22
ln -s ${InputDir_Liege}/220816_A00801_0259_BHJGTNDSX3_20220812-GEN0911-PL23-24/fastq_bulk_All ${OutputDir}/plate-23
ln -s ${InputDir_Liege}/220816_A00801_0259_BHJGTNDSX3_20220812-GEN0911-PL25-26/fastq_bulk_All ${OutputDir}/plate-25
ln -s ${InputDir_Liege}/221011_A00801_0278_BH3MM2DSX5_20220921-Hperee-GEN1033-P27-28/fastq_bulk_All ${OutputDir}/plate-27
ln -s ${InputDir_Liege}/221011_A00801_0278_BH3MM2DSX5_20220921-Hperee-GEN1033-P29-30/fastq_bulk_All ${OutputDir}/plate-29
ln -s ${InputDir_Liege}/221011_A00801_0278_BH3MM2DSX5_20220921-Hperee-GEN1033-P31-32/fastq_bulk_All ${OutputDir}/plate-31
ln -s ${InputDir_Liege}/230203_A00801_0332_AHCGWLDSX5_20230130-GEN1193-Hperee-PL33/fastq_bulk_All ${OutputDir}/plate-33
ln -s ${InputDir_Liege}/230203_A00801_0332_AHCGWLDSX5_20230130-GEN1193-Hperee-PL34/fastq_bulk_All ${OutputDir}/plate-34
ln -s ${InputDir_Liege}/230419_A00801_0364_AHNG5KDSX5_20230417-GEN1273-Hperee-Pl35-36/fastq_bulk_All ${OutputDir}/plate-35
ln -s ${InputDir_Liege}/230203_A00801_0332_AHCGWLDSX5_20230130-GEN1193-Hperee-PL37/fastq_bulk_All ${OutputDir}/plate-37
ln -s ${InputDir_Liege}/230203_A00801_0332_AHCGWLDSX5_20230130-GEN1193-Hperee-PL38/fastq_bulk_All ${OutputDir}/plate-38
ln -s ${InputDir_Liege}/230419_A00801_0364_AHNG5KDSX5_20230417-GEN1273-Hperee-Pl39/fastq_bulk_All ${OutputDir}/plate-39
ln -s ${InputDir_Liege}/230419_A00801_0364_AHNG5KDSX5_20230417-GEN1273-Hperee-Pl40/fastq_bulk_All ${OutputDir}/plate-40
ln -s ${InputDir_Liege}/230427_A00801_0366_AHVNKTDSX5_20230425-GEN1338-Hperee-Plate41-42/fastq_bulk_All ${OutputDir}/plate-41
ln -s ${InputDir_Liege}/230427_A00801_0366_AHVNKTDSX5_20230425-GEN1338-Hperee-Plate43-44/fastq_bulk_All ${OutputDir}/plate-43
ln -s ${InputDir_Liege}/230427_A00801_0366_AHVNKTDSX5_20230425-GEN1338-Hperee-Plate45-46/fastq_bulk_All ${OutputDir}/plate-45
ln -s ${InputDir_Liege}/230602_A00801_0386_AH27KHDSX7_20230530-GEN1375-Peree-Plate47-48/fastq_bulk_All ${OutputDir}/plate-47
ln -s ${InputDir_Liege}/230602_A00801_0386_AH27KHDSX7_20230530-GEN1375-Peree-Plate49-50 ${OutputDir}/plate-49
ln -s ${InputDir_Liege}/230623_A00801_0393_AH2NWMDSX7_20230620-Hperee-plate51/fastq_bulk_All ${OutputDir}/plate-51
ln -s ${InputDir_Liege}/230623_A00801_0393_AH2NWMDSX7_20230620-GEN1377-Hperee-Plate52 ${OutputDir}/plate-52
ln -s ${InputDir_Liege}/230623_A00801_0393_AH2NWMDSX7_20230620-GEN1377-Hperee-Plate53 ${OutputDir}/plate-53
ln -s ${InputDir_Liege}/230623_A00801_0393_AH2NWMDSX7_20230620-GEN1377-Hperee-Plate54/fastq_bulk_All ${OutputDir}/plate-54
ln -s ${InputDir_Liege}/230825_A00801_0420_BHHC3MDSX7_20230825-GEN1409-Hperee-PL55/fastq_bulk_L1-3 ${OutputDir}/plate-55
ln -s ${InputDir_Liege}/230825_A00801_0420_BHHC3MDSX7_20230825-GEN1409-Hperee-Pl56/fastq_bulk_L1-3 ${OutputDir}/plate-56
ln -s ${InputDir_Liege}/230825_A00801_0420_BHHC3MDSX7_20230825-GEN1409-Hperee-Pl57/fastq_bulk_L1-3 ${OutputDir}/plate-57
ln -s ${InputDir_Liege}/230825_A00801_0420_BHHC3MDSX7_20230825-GEN1409-Hperee-Pl58 ${OutputDir}/plate-58
ln -s ${InputDir_Liege}/230915_A00801_0433_BHKYVKDSX7_20230908-GEN1482-Hperee-Plate59/fastq_Bulk ${OutputDir}/plate-59
ln -s ${InputDir_Liege}/230915_A00801_0433_BHKYVKDSX7_20230908-GEN1482-Hperee-Plate60/fastq_Bulk ${OutputDir}/plate-60
ln -s ${InputDir_Liege}/230915_A00801_0433_BHKYVKDSX7_20230907-GEN1482-Hperee-plate61/fastq_Bulk ${OutputDir}/plate-61
ln -s ${InputDir_Liege}/230915_A00801_0433_BHKYVKDSX7_20230908-GEN1482-Hperee-PLate62/fastq_Bulk ${OutputDir}/plate-62
ln -s ${InputDir_Liege}/231120_A00801_0465_AHJ52YDSX7_20231120-GEN1504-Hperee-P63/fastq_bulk_All ${OutputDir}/plate-63
ln -s ${InputDir_Liege}/231120_A00801_0465_AHJ52YDSX7_20231120-GEN1504-Hperee-P64/fastq_bulk_All ${OutputDir}/plate-64
ln -s ${InputDir_Liege}/240206_A00801_0501_AH5G5HDSXC_20240201-GEN1581-Hperee-P65/fastq_Bulk ${OutputDir}/plate-65
ln -s ${InputDir_Liege}/240206_A00801_0501_AH5G5HDSXC_20240201-GEN1581-Hperee-P66/fastq_Bulk ${OutputDir}/plate-66
ln -s ${InputDir_Liege}/231120_A00801_0465_AHJ52YDSX7_20231120-GEN1504-Hperee-P67/fastq_bulk_All ${OutputDir}/plate-67
ln -s ${InputDir_Liege}/231120_A00801_0465_AHJ52YDSX7_20231120-GEN1504-Hperee-P68/fastq_bulk_All ${OutputDir}/plate-68
ln -s ${InputDir_Liege}/240802_A00801_0594_AHF7FWDSXC_20240716-GEN1916-Hperee-Plate69-70/fastq_Bulk ${OutputDir}/plate-69
ln -s ${InputDir_Liege}/240508_A00801_0553_AHCY7CDSXC_20240430-Hperee-P71/fastq_Bulk ${OutputDir}/plate-71
ln -s ${InputDir_Liege}/240508_A00801_0553_AHCY7CDSXC_20240430-Hperee-P72/fastq_Bulk ${OutputDir}/plate-72
ln -s ${InputDir_Liege}/240412_A00937_0079_AHHWCNDSX7_20240411-GEN1707-Hperee-Pl73/fastq_Bulk ${OutputDir}/plate-73
ln -s ${InputDir_Liege}/240412_A00937_0079_AHHWCNDSX7_20240411-GEN1707-Hperee-Pl74/fastq_Bulk ${OutputDir}/plate-74
ln -s ${InputDir_Liege}/240206_A00801_0501_AH5G5HDSXC_20240201-GEN1581-Hperee-P75 ${OutputDir}/plate-75
ln -s ${InputDir_Liege}/240206_A00801_0501_AH5G5HDSXC_20240201-GEN1581-Hperee-P76 ${OutputDir}/plate-76
ln -s ${InputDir_Liege}/240412_A00937_0079_AHHWCNDSX7_20240411-GEN1707-Hperee-Pl77/fastq_Bulk ${OutputDir}/plate-77
ln -s ${InputDir_Liege}/240412_A00937_0079_AHHWCNDSX7_20240411-GEN1707-Hperee-Pl78/fastq_Bulk ${OutputDir}/plate-78
ln -s ${InputDir_Liege}/240606_A00801_0567_AHF77LDSXC_20240430-Hperee-P79/fastq_Bulk ${OutputDir}/plate-79
ln -s ${InputDir_Liege}/240606_A00801_0567_AHF77LDSXC_20240430-Hperee-P80/fastq_Bulk ${OutputDir}/plate-80
ln -s ${InputDir_Liege}/240802_A00801_0594_AHF7FWDSXC_20240716-GEN1916-Hperee-Plate81-82/fastq_Bulk ${OutputDir}/plate-81
ln -s ${InputDir_Liege}/240628_A00801_0579_AHHW2JDSXC_20240627-GEN1890-P83-84/fastq_Bulk ${OutputDir}/plate-83
ln -s ${InputDir_Liege}/240628_A00801_0579_AHHW2JDSXC_20240627-GEN1890-P85-86/fastq_Bulk ${OutputDir}/plate-85
ln -s ${InputDir_Liege}/240802_A00801_0594_AHF7FWDSXC_20240716-GEN1916-Hperee-plate87-88/fastq_Bulk ${OutputDir}/plate-87

echo -e "\n"
echo "Symlinking files completed"
echo -e "\n"
date
