library(Rsamtools)
library(readr)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))



## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data","independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load snv_df file =====================
# OK tell me what kind of SNV mutation you like
# you can choose all or driver
# write it down
snvtype <-  "driver" # you can write down driver or all
if (snvtype == "driver") {
  snv_df<-read.csv(file.path(root_dir,"scratch","sv-snv","driver_snv_in_samples_and_genes.csv"),row.names=1)
} else {
  snv_df<-read.csv(file.path(root_dir,"scratch","sv-snv","snv_in_samples_and_genes.csv"),row.names=1)
}


## ===================== Load breakpoint =========================
breakpoint <- read_tsv(file.path(root_dir,"data","pbta-sv-manta.tsv.gz"))
breakpoint_bnd <- breakpoint[breakpoint$SV.type == "BND"& breakpoint$FILTER == "PASS",c("SV.chrom","SV.start","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")]
breakpoint_unbnd_start <- breakpoint[breakpoint$SV.type != "BND"& breakpoint$FILTER == "PASS",c("SV.chrom","SV.start","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")]
breakpoint_unbnd_end <- breakpoint[breakpoint$SV.type != "BND"& breakpoint$FILTER == "PASS",c("SV.chrom","SV.end","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")]
colnames(breakpoint_bnd) <- c("SV.chrom","breakpoint","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")
colnames(breakpoint_unbnd_start) <- c("SV.chrom","breakpoint","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")
colnames(breakpoint_unbnd_end) <- c("SV.chrom","breakpoint","SV.length","SV.type","Kids.First.Biospecimen.ID.Tumor")
breakpoint_all <- rbind(breakpoint_bnd,
                        breakpoint_unbnd_start,
                        breakpoint_unbnd_end)

## ==================== Divided patients into 2 groups ============
mygene <- "H3F3A"
mutation_bioid <- colnames(snv_df[,snv_df[mygene,] > 0])
nonmutation_bioid <-  bioid[(bioid %in% mutation_bioid) == FALSE]


## ==================== Count number and do test=================
# FUNCTION find in and out site
in_out <- function(samplesv)  {
  conclusion_df <- data.frame()
  for (i in 1:nrow(samplesv)){
    conclusion = sum((bedfile$X1 == samplesv[i,1]) & (samplesv[i,2] >= bedfile$X2) &  (samplesv[i,2] <= bedfile$X3))
    conclusion_df = rbind(conclusion_df,conclusion)
  }
  colnames(conclusion_df) = "chipseq"
  return(conclusion_df)
}


# read chipseq bed file

allfiles <- list.files(file.path(root_dir,"scratch","chipseq"),pattern = "*.gz")
resultout <- data.frame(matrix(NA,ncol = 7,nrow = 0))
colnames(resultout) <- c("breakpoint_in_mutaion",
                         "breakpoint_in_unmutaion",
                         "breakpoint_all_mutaion",
                         "breakpoint_all_unmutaion",
                         "p.value",
                         "or",
                         "file"
                         )
for (filei in allfiles) {
  
  filepath <- file.path(root_dir,"scratch", "chipseq", filei)
  bedfile <- read_tsv(filepath, col_names=F)
  bedfile$X1 <- gsub("chr","",bedfile$X1)
  
  breakpoint_inout_mutaion <- data.frame()
  for (i  in mutation_bioid) {
    samplesv <- as.data.frame(breakpoint_all[breakpoint_all$Kids.First.Biospecimen.ID.Tumor == i,])
    samplesv = cbind(samplesv,in_out(samplesv))
    breakpoint_inout_mutaion = rbind(breakpoint_inout_mutaion,samplesv)
  }
  breakpoint_inout_unmutaion <- data.frame()
  for (i  in nonmutation_bioid) {
    samplesv <- as.data.frame(breakpoint_all[breakpoint_all$Kids.First.Biospecimen.ID.Tumor == i,])
    samplesv = cbind(samplesv,in_out(samplesv))
    breakpoint_inout_unmutaion = rbind(breakpoint_inout_unmutaion,samplesv)
  }
  nrow(breakpoint_inout_mutaion[breakpoint_inout_mutaion$chipseq>0,])/nrow(breakpoint_inout_mutaion)
  nrow(breakpoint_inout_unmutaion[breakpoint_inout_unmutaion$chipseq>0,])/nrow(breakpoint_inout_unmutaion)
  result <- fisher.test(matrix(data=c(nrow(breakpoint_inout_mutaion[breakpoint_inout_mutaion$chipseq>0,]),
                                      nrow(breakpoint_inout_mutaion),
                                      nrow(breakpoint_inout_unmutaion[breakpoint_inout_unmutaion$chipseq>0,]),
                                      nrow(breakpoint_inout_unmutaion)),nrow = 2,ncol = 2))
  if (result$p.value < 0.05 & result$estimate > 1) {
    resultadd <- data.frame("breakpoint_in_mutaion"=nrow(breakpoint_inout_mutaion[breakpoint_inout_mutaion$chipseq>0,]),
                            "breakpoint_in_unmutaion"=nrow(breakpoint_inout_unmutaion[breakpoint_inout_unmutaion$chipseq>0,]),
                            "breakpoint_all_mutaion"=nrow(breakpoint_inout_mutaion),
                            "breakpoint_all_unmutaion"=nrow(breakpoint_inout_unmutaion),
                            "p.value"=result$p.value,
                            "or"=result$estimate,
                            "file"=filei)
    resultout <- rbind(resultout,resultadd)
  }
}
write.csv(resultout,file.path(root_dir,"scratch","chipseq","chipseq_significant_result.csv"),quote = F)


## ===================== Load chipseq info ========================
chipseq_info <- read.delim(file.path(root_dir,"scratch","chipseq","experiment_report_2020_2_10_17h_44m.tsv"))
resultout_name  <- gsub(".bed.gz","",unique(resultout$file))
allfiles_name <- gsub(".bed.gz","",list.files(file.path(root_dir,"scratch","chipseq"),pattern = "*.gz"))

resultout_info <- data.frame(matrix(NA,ncol =8 ,nrow = 0))    
colnames(resultout_info) <-  c("name","Target.of.assay","Biosample.summary","Biosample.term.name","Life.stage","Age","Age.units","Biosample.treatment")
for (namei in resultout_name){
  resultout_add <- chipseq_info[grepl(namei, chipseq_info$Files),c("Target.of.assay","Biosample.summary","Biosample.term.name","Life.stage","Age","Age.units","Biosample.treatment")]
  resultout_add$name  <- namei
  resultout_info <- rbind(resultout_info,resultout_add)
}
allfiles_info <- data.frame(matrix(NA,ncol =8 ,nrow = 0))    
colnames(allfiles_info) <-  c("name","Target.of.assay","Biosample.summary","Biosample.term.name","Life.stage","Age","Age.units","Biosample.treatment")
for (namei in allfiles_name){
  resultout_add <- chipseq_info[grepl(namei, chipseq_info$Files),c("Target.of.assay","Biosample.summary","Biosample.term.name","Life.stage","Age","Age.units","Biosample.treatment")]
  resultout_add$name  <- namei
  allfiles_info <- rbind(allfiles_info,resultout_add)
}



# combine resultout_info and resultout
resultout$name <- gsub(".bed.gz","",resultout$file)
resultout_info <- merge(resultout,resultout_info)

# WRITE
write.csv(resultout_info,file.path(root_dir,"scratch","chipseq","H3 chipseq significant result and sample info.csv"),quote = F)
