# Usage example - from Inf_Fibrosarcoma project 
setwd("/mnt/scratch1/maycon/ToolBox/BulkRNAseq_Processing_STAR")

# On terminal --------
# First, Activate the conda environment which has STAR all set up
# Run the following command: Rscript /mnt/scratch1/maycon/ToolBox/BulkRNAseq_Processing_STAR/Alignment.R & 



# The script running on background ---------
# 1. Define variables  
files <- dir("/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/", pattern = "fastq")
path = "/media/ResearchHome/plummgrp/home/common/LFurtado-colab/RNA/FASTQ_files_ETV6_NTRK3_fused_tumors/"
path.out = "/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/"

unique_sample_ids <- c("SJBT032721_D1-S3", 
                       "SJIFS031085_D1-S9",
                       "SJST030433_D1-S3", 
                       "SJST030567_D1-S10", 
                       "SJST032838_D1-S4",  
                       "SJST032952_D1-S9", 
                       "SJST032952_D2-S5", 
                       "SJST032952_D4-S4", 
                       "SJST032952_D4-S9",  
                       "SJST033312_D1-S13", 
                       "SJST033312_D1-S2") 

# 2. Running STAR pipeline 
# library(dplyr)
for(i in 1:length(unique_sample_ids)) {
  #i <- 1
  samples <- files[grep(unique_sample_ids[i], files)]
  samples_R1 <- samples[grep("R1_0", samples)]
  samples_R2 <- samples[grep("R2_0", samples)]
  
  if (length(samples) > 2) { #samples with multiple runs 
    system(paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
                  path, samples_R1[1],',',path,samples_R1[2],' ',path,samples_R2[1],',',path,samples_R2[2], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i]))
    
  } else { #usual 2 runs samples (2 fastq)
    system((paste0("STAR --genomeDir /media/ResearchHome/plummgrp/home/common/GRCh37_homosapiens/HS_STAR_dir_2 --runThreadN 20 --readFilesIn ",
                   path, samples_R1[1],' ',path,samples_R2[1], " --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path.out, unique_sample_ids[i])))
  }
}


# # On Rstudio ------------
# # Check it if we ran all fastq samples 
# # all the fastq samples
# unique_sample_ids
# # catch all .bam files 
# library(stringr)
# bam_out <- dir("/mnt/scratch1/LFurtado-fibrosarcoma/RNA-Align/ETV6_smps/", pattern = ".bam$")
# file_split_1 <- str_split_fixed(as.character(bam_out), "[Aligned]", 2)[,1]
# unique_sample_ids_check <- unique(file_split_1)
# 
# # If it's TRUE so you're all covered. If it's FALSE, take a look which of the samples didn't go through 
# identical(unique_sample_ids, unique_sample_ids_check) #FALSE 
# # SJST032952_D1-S9 - NEED TO REVIEW THE FILE PAIR
# # SJST032952_D4-S9 - NEED TO REVIEW - ONLY ONE FASTQ FILE. It has a .bam output but it's not right
# # All the rest went through! 


