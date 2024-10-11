#!/bin/bash

#SBATCH --partition=kosmos
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --job-name=sc_eqtl_1
#SBATCH --output=logs/nf_dsl2_sc_eqtl_%j.log

module load jdk/jdk1.8.0_66
module load singularity/3.2.1
module load nextflow/21.03.0-edge
module load slurm


SingularityDir=.


export _JAVA_OPTIONS=-Djava.io.tmpdir=.
export NXF_SINGULARITY_CACHEDIR=${SingularityDir}

srun nextflow run nextflow_eQTL_dsl2_1.nf -c nextflow.config -resume -with-timeline pipeline_timeline.html -with-trace --number_of_peer_iterations 2000 --peer_approach 


echo "" 
echo "############################################################################"
echo "sacct output:" 
echo "" 
sacct --format="JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible" -j ${SLURM_JOB_ID} 
echo ""
date
echo "zhmi zaletnye!"


