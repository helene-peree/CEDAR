

data_folder_cluster <- ""
data_folder <- ""

sequencing_metadata <- read.csv(paste0(data_folder,"metadata_sequencing/50_samples.csv"))

sequencing_metadata$fastq.file.cDNA <- sapply(sequencing_metadata$fastq.file.cDNA, function(x){
  
  sub("~","/home/u/",x, fixed = T)
  
})

sequencing_metadata$cDNA.run.name <- sapply(sequencing_metadata$cDNA.run.name, function(x){
  
  sub("_S.*","",x, fixed = F)
  
})

sequencing_metadata$fastq.file.HTO. <- sapply(sequencing_metadata$fastq.file.HTO., function(x){
  
  sub("~","/home/u/",x, fixed = T)
  
})

sequencing_metadata$HTO.run.name <- sapply(sequencing_metadata$HTO.run.name, function(x){
  
  sub("_S.*","",x, fixed = F)
  
})

for (i in 1:nrow(sequencing_metadata)){
  
  libraries <- as.data.frame(matrix(ncol = 3, nrow = 2))
  colnames(libraries) <- c("fastqs","sample","library_type")
  
  libraries[1,1] <- sequencing_metadata$fastq.file.cDNA[i]
  libraries[1,2] <- sequencing_metadata$cDNA.run.name[i]
  libraries[2,1] <- sequencing_metadata$fastq.file.HTO.[i]
  libraries[2,2] <- sequencing_metadata$HTO.run.name[i]
  libraries[,3] <- c("Gene Expression","Antibody Capture")
  
  write.csv(libraries,paste0(data_folder,"metadata_sequencing/libraries/",
                             sequencing_metadata$sample[i],"_lib.csv"), quote = F, row.names = F)
  
  #####################_make_script_############################
  
  script_name <- paste0(data_folder, "scripts/sample_", 
                sequencing_metadata$sample[i],"_script.sh")
  
  file.create(script_name)
  
  sink(script_name)
  cat("#!/bin/bash")
  cat("\n")
  cat("\n")
  cat("#SBATCH --partition=all_24hrs,kosmos")
  cat("\n")
  cat("#SBATCH --ntasks=1")
  cat("\n")
  cat("#SBATCH --cpus-per-task=10")
  cat("\n")
  cat("#SBATCH --mem-per-cpu=14G")
  cat("\n")
  cat(paste0("#SBATCH --job-name=cellranger_", sequencing_metadata$sample[i]))
  cat("\n")
  cat(paste0("#SBATCH --output=",data_folder_cluster,
             "logs/cellr_",sequencing_metadata$sample[i],"_%j.log"))
  cat("\n")
  cat("\n")
  cat("module load cellranger/7.1.0")
  cat("\n")
  cat("\n")
  cat("cd ../cellranger_results")
  cat("\n")
  cat(paste0("mkdir sample_",sequencing_metadata$sample[i]))
  cat("\n")
  cat(paste0("cd sample_",sequencing_metadata$sample[i]))
  cat("\n")
  cat("\n")
  cat(paste0("cellranger count --id=",sequencing_metadata$sample[i],
             " --transcriptome=/Resources/Genomes/CellRanger/Homo_sapiens/Ensembl/GRCh38/release_103 ",
             "--libraries=", data_folder_cluster,"metadata_sequencing/libraries/",
             sequencing_metadata$sample[i],"_lib.csv --feature-ref=",
             data_folder_cluster,"metadata_sequencing/TotalSeqB_feature_ref.csv"))
  cat("\n")
  cat("\n")
  cat("sacct --format='JobId,JobName,NodeList,State,Elapsed,Timelimit,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveVMSize,ReqMem,Submit,Eligible' -j ${SLURM_JOB_ID}")
  cat("\n")
  cat("\n")
  
  sink()
  
}

file.create(paste0(data_folder, "scripts/ruler_sc.sh"))

sink(paste0(data_folder, "scripts/ruler_sc.sh"))

cat("#!/bin/bash")
cat("\n")
cat("\n")
cat("#SBATCH --partition=all_24hrs,all_5hrs")
cat("\n")
cat("#SBATCH --ntasks=1")
cat("\n")
cat("#SBATCH --mem-per-cpu=4G")
cat("\n")
cat(paste0("#SBATCH --job-name=ruler-%A_%a.out"))
cat("\n")
cat(paste0("#SBATCH --output=",data_folder_cluster,
           "logs/rul_out_%j.log"))
cat("\n")
cat("\n")
cat("module load slurm")
cat("\n")
cat("\n")
cat(paste0("job_1=$(sbatch sample_B11_script.sh | awk '{print $4}')"))
cat("\n")
for (i in 2:nrow(sequencing_metadata)){
  
  cat(paste0("job_",i,
             '=$(sbatch --dependency=afterok:"$job_',i-1,'" sample_',
             sequencing_metadata$sample[i],"_script.sh | awk '{print $4}')"))
  cat("\n")
}

cat("\n")

sink()














