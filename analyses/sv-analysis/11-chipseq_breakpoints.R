library(Rsamtools)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# read chipseq bed file
filepath <- file.path(root_dir,"scratch", "chipseq", "ENCFF436ITH.bed.gz")
bedfile <- read_tsv(filepath, col_names=F)
bedfile$X1 <- gsub("chr","",bedfile$X1)

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


## ==================== Count number =================
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

breakpoint_inout_mutaion <- data.frame()
for (i  in mutation_bioid) {
  print(i)
  samplesv <- as.data.frame(breakpoint_all[breakpoint_all$Kids.First.Biospecimen.ID.Tumor == i,])
  samplesv = cbind(samplesv,in_out(samplesv))
  breakpoint_inout_mutaion = rbind(breakpoint_inout_mutaion,samplesv)
}
breakpoint_inout_unmutaion <- data.frame()
for (i  in nonmutation_bioid) {
  print(i)
  samplesv <- as.data.frame(breakpoint_all[breakpoint_all$Kids.First.Biospecimen.ID.Tumor == i,])
  samplesv = cbind(samplesv,in_out(samplesv))
  breakpoint_inout_unmutaion = rbind(breakpoint_inout_unmutaion,samplesv)
}
nrow(breakpoint_inout_mutaion[breakpoint_inout_mutaion$chipseq>0,])/nrow(breakpoint_inout_mutaion)
nrow(breakpoint_inout_unmutaion[breakpoint_inout_unmutaion$chipseq>0,])/nrow(breakpoint_inout_unmutaion)
fisher.test(matrix(data=c(nrow(breakpoint_inout_mutaion[breakpoint_inout_mutaion$chipseq>0,]),
                          nrow(breakpoint_inout_mutaion),
                          nrow(breakpoint_inout_unmutaion[breakpoint_inout_unmutaion$chipseq>0,]),
                          nrow(breakpoint_inout_unmutaion)),nrow = 2,ncol = 2))
