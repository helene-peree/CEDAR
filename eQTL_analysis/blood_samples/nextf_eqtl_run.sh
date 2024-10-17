#!/bin/bash
#
#SBATCH --partition=nextflow
#
#SBATCH --job-name=nf_eQTL
#SBATCH --output=../slurm_logs/nf_eQTL_%j.log
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36G

module load EasyBuild
module load singularity/3.7.1
module load Nextflow/22.10.3
module load R/4.2.2-foss-2022b
module load java/8.0.111
module load BCFtools/1.17-GCC-12.2.0

ScratchDir="/scratch/GIGA/USER/u230859"
WorkingDir=$1
nb_factors=$2

cd ${WorkingDir}
export _JAVA_OPTIONS=-Djava.io.tmpdir=${ScratchDir}
export NXF_SINGULARITY_CACHEDIR=${WorkingDir}

nextflow run nextflow_eQTL_dsl2_beta1.2.2.nf -c nextflow.config --outdir ${ScratchDir} -resume -with-timeline pipeline_timeline.html -with-trace --No_peer_factors ${nb_factors}

echo "" 
echo "############################################################################"
echo "sacct output:" 
echo "" 
sacct --format="JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible" -j ${SLURM_JOB_ID} 
echo ""
date
echo "Cis eQTL mapping completed"
