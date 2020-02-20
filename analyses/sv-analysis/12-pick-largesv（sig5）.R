if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.masked")


setwd("~/GitHub/OpenPBTA-analysis")

library(ConsensusClusterPlus)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(ggplot2)
library(ComplexHeatmap)
library(readr)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_directory <- file.path(root_dir, "scratch","sv-largeall")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data","independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Read simple SV and Keep Large Del/Dup/Inv, cluster_H3mutation_Disease===================
# have done
# can skip
largeallfun <- function(simple_sv){
  size <-  abs(as.numeric(simple_sv["pos1"])-as.numeric(simple_sv["pos2"]))
  out <- ifelse(simple_sv[7] %in% c("DEL","DUP","t2tINV","h2hINV") & (size > 1000000),"largeall","nolargeall")
  return(out)
}

largeall_stat <- data.frame(matrix(NA,nrow = 0,ncol = 3))
colnames(largeall_stat) <- c("bioid","SimpleSVNum","LargeAllNum")
for (bioidi in bioid) {
  simple_sv <- read.delim(file.path(root_dir,"scratch","sv-shatterseek","tsv",paste0(bioidi,"_removechrss.tsv")))
  simple_sv$largeall <- apply(simple_sv,1,largeallfun)
  write.table(simple_sv,file.path(output_directory,paste0(bioidi,"_largeall.tsv")),sep = "\t",row.names = F,col.names = T)
  largeall_stat_add <- data.frame("bioid"=bioidi,
                                  "SimpleSVNum"=nrow(simple_sv),
                                  "LargeAllNum"=nrow(simple_sv[simple_sv$largeall=="largeall",]))
  largeall_stat <- rbind(largeall_stat,largeall_stat_add)
}


# load SV signatures
sv_signature_path <- file.path("D:","Yang","Downloads","SigProfiler_2_5_1_7","SigProfiler_2_5_1_7","output","PBTA_PASS_v7.112","text","res_PBTA_PASS_v7.112_signature_activities_for_5_sigs.csv")
sv_signature <- read.csv(sv_signature_path)
sv_signature$Sample.Names <- paste0("BS_",sv_signature$Sample.Names)
row.names(sv_signature) <- sv_signature$Sample.Names
sv_signature <-t(sv_signature[,-1,drop=FALSE])
# Cluster:consensusClass
title <- "SV Signatures Clusters"
# do cluster
results <- ConsensusClusterPlus(sv_signature,maxK = 7,reps = 50,pItem=0.8,pFeature = 1,
                                title = title, clusterAlg = "hc",distance = "pearson", seed = 123456,plot = "png")

#  choose cluster number
clusternumber <- 5
#  add cluster result to sv_signature
sv_signature <- rbind(sv_signature,results[[clusternumber]][["consensusClass"]])
row.names(sv_signature)[nrow(sv_signature)] <- "consensusClass"
sv_signature <- as.data.frame(t(sv_signature))
# combine cluster
sv_signature$bioid <- row.names(sv_signature)
largeall_stat <- merge(largeall_stat,sv_signature[,c("bioid","consensusClass")])


# Load H3 mutation
# read snv file
snv <- read.csv(file.path(root_dir,"scratch","sv-snv","driver_snv_in_samples_and_genes.csv"))
# input mygene
mygene <- "H3F3A"
mygene_snv <- snv[snv$X == mygene,]
mygene_snv <- as.data.frame(t(mygene_snv[,2:ncol(mygene_snv)]))
colnames(mygene_snv) <- mygene
mygene_snv$bioid <- row.names(mygene_snv)
# combine mygene mutation
largeall_stat <- merge(largeall_stat,mygene_snv,all.x = T)
largeall_stat[is.na(largeall_stat)] <- 0

# Load Disease Info
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
# combine disease
largeall_stat <- merge(largeall_stat,info[,c("Kids_First_Biospecimen_ID","short_histology")],by.x = "bioid",by.y = "Kids_First_Biospecimen_ID",all.x = T)

# write out
write.csv(largeall_stat,file.path(output_directory,"patients_largeall_cluster_H3mutation_disease.csv"),quote = F,row.names = F)

## ====================== Load Gap =================
# no SV fall into heterochromatin and telomere
# Skip this section

# function ingap
ingap <- function(simple_sv,gapfile){
  chr1 <- paste0("chr",simple_sv["chrom1"])
  pos1 <- as.numeric(simple_sv["pos1"])
  chr2 <- paste0("chr",simple_sv["chrom2"])
  pos2 <- as.numeric(simple_sv["pos2"])
  if (sum((chr1 == gapfile$chrom) & (pos1 >= as.numeric(gapfile$chromStart)) & (pos1 <= as.numeric(gapfile$chromEnd)))>0) {
    out <-  "yes"
    print("ha")
  } else {
    out <- "na"
  }
  if(sum((chr2 == gapfile$chrom) & (pos2 >= as.numeric(gapfile$chromStart)) & (pos2 <= as.numeric(gapfile$chromEnd)))>0) {
    out <-  paste(out,"yes",sep = "_")
    print("ha")
  } else {
    out <- paste(out,"na",sep = "_")
  }
  return(out)
  
}


gap <- read.delim(file.path(root_dir,"scratch","chipseq","hg38_gap"))
gap_telomere <- gap[gap$type == "telomere",]
gap_heterochromatin <- gap[gap$type == "heterochromatin",]
# add in gap
for (bioidi in bioid) {
  simple_sv <- read.delim(file.path(root_dir,"scratch","sv-largeall",paste0(bioidi,"_largeall.tsv")))
  simple_sv$telomere <- apply(simple_sv, 1, ingap,gapfile=gap_telomere)
  simple_sv$heterochromatin <- apply(simple_sv, 1, ingap,gapfile=gap_heterochromatin)
  write.table(simple_sv,file.path(output_directory,paste0(bioidi,"_largeall.tsv")),sep = "\t",row.names = F,col.names = T)
}


## ===================== Plot Large SV breakpoint ===================
genome <- BSgenome.Hsapiens.UCSC.hg38.masked
length <- c(length(genome$chr1),length(genome$chr2),length(genome$chr3),
            length(genome$chr4),length(genome$chr5),length(genome$chr6),
            length(genome$chr7),length(genome$chr8),length(genome$chr9),
            length(genome$chr10),length(genome$chr11),length(genome$chr12),
            length(genome$chr13),length(genome$chr14),length(genome$chr15),
            length(genome$chr16),length(genome$chr17),length(genome$chr18),
            length(genome$chr19),length(genome$chr20),length(genome$chr21),
            length(genome$chr22),length(genome$chrX),length(genome$chrY))
length_sum <- c()
for  (i in 1:24) {
  length_sum <- c(length_sum,sum(length[1:i]))
}
length_sum2 <- c()
for  (i in 1:24) {
  length_sum2 <- c(length_sum2,sum(length[0:(i-1)]))
}


# function record location
breakpointlocation <- function(simple_sv,pos){
  chrnumber <- simple_sv["chrom1"]
  if (chrnumber == "X") {
    chrnumber <- 23
  }
  if (chrnumber == "Y") {
    chrnumber <- 24
  }
  position <- as.numeric(simple_sv[pos]) 
  location <- sum(length[1:chrnumber-1])+position
  return(location)
}
# function get chr
locationchr <- function(coordinate) {
  out <- min(which( (length_sum - coordinate[1])>=0 ))
}


# record all breakpoint's location
largeall_breakpoint <- c()
for (bioidi in bioid) {
  simple_sv <- read.delim(file.path(root_dir,"scratch","sv-largeall",paste0(bioidi,"_largeall.tsv")))
  simple_sv <- simple_sv[simple_sv$largeall == "largeall",]
  if (nrow(simple_sv) > 0){
    breakpoint <- as.numeric(apply(simple_sv, 1, breakpointlocation,pos="pos1"))
    largeall_breakpoint <- c(largeall_breakpoint,breakpoint)
    breakpoint <- as.numeric(apply(simple_sv, 1, breakpointlocation,pos="pos2"))
    largeall_breakpoint <- c(largeall_breakpoint,breakpoint)
  }
}
# record nonlargeall breakpoint's location
nonlargeall_breakpoint <- c()
for (bioidi in bioid) {
  simple_sv <- read.delim(file.path(root_dir,"scratch","sv-largeall",paste0(bioidi,"_largeall.tsv")))
  simple_sv <- simple_sv[simple_sv$largeall != "largeall" & simple_sv$chrom1 != "M",]
  if (nrow(simple_sv) > 0){
    breakpoint <- as.numeric(apply(simple_sv, 1, breakpointlocation,pos="pos1"))
    nonlargeall_breakpoint <- c(nonlargeall_breakpoint,breakpoint)
    breakpoint <- as.numeric(apply(simple_sv, 1, breakpointlocation,pos="pos2"))
    nonlargeall_breakpoint <- c(nonlargeall_breakpoint,breakpoint)
  }
}



# build a genome coordinate
coordinate <- data.frame(matrix(seq(0,sum(length[1:24]),1000000)))
a=table(cut(largeall_breakpoint, breaks = seq(0,sum(length[1:24]),1000000)))
coordinate$freq <- c(a,0)
coordinate$chrom <- apply(coordinate,1,locationchr)
coordinate$chrom <- as.factor(coordinate$chrom)
# all nonlargeall breakpoint freq
a=table(cut(nonlargeall_breakpoint, breaks = seq(0,sum(length[1:24]),1000000)))
coordinate$freq_nonlargeall <- c(a,0)
# add position
coordinate$position <- as.numeric(coordinate[,1]) - length_sum2[coordinate$chrom]


# plot
lower <- HeatmapAnnotation("chrom" = coordinate[,"chrom"],
                           show_legend = F)

upper <- HeatmapAnnotation(
  freqnumber_largeall = anno_barplot(coordinate[,"freq"]),
  freqnumber_nonlargeall = anno_barplot(coordinate[,"freq_nonlargeall"]),
  height = unit(5, "cm"),
  width = 1,
  show_legend = TRUE
)

zero_row_mat <- matrix(nrow = 0, ncol = nrow(coordinate))

ht <- Heatmap(zero_row_mat,
              top_annotation = upper,
              bottom_annotation = lower,
              show_heatmap_legend = FALSE,
              show_column_names = FALSE)
draw(ht)



## =================== Load chipseq files ==================
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

# summarize all sv
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
for ( i in  1:nrow(breakpoint_all)) {
  if (breakpoint_all[i,"SV.type"] != "BND"  & breakpoint_all[i,"SV.length"] > 1000000) {
    breakpoint_all[i,"largeall"] <-  "largeall"
  }  else {
    breakpoint_all[i,"largeall"] <-  "nolargeall"
  }
}


# read chipseq bed file
allfiles <- list.files(file.path(root_dir,"scratch","chipseq"),pattern = "*.gz")
resultout <- data.frame(matrix(NA,ncol = 7,nrow = 0))
colnames(resultout) <- c("largeall_in",
                         "largeall_out",
                         "nolargeall_in",
                         "nolargeall_out",
                         "p.value",
                         "or",
                         "file"
)

count <- 0
for (filei in allfiles) {
  count <- count +1
  cat(count)
  
  filepath <- file.path(root_dir,"scratch", "chipseq", filei)
  bedfile <- read_tsv(filepath, col_names=F)
  bedfile$X1 <- gsub("chr","",bedfile$X1)
  

  
  samplesv_largeall <- as.data.frame(breakpoint_all[breakpoint_all$largeall == "largeall",])
  samplesv_largeall = cbind(samplesv_largeall,in_out(samplesv_largeall))
  
  
  samplesv_nolargeall <- as.data.frame(breakpoint_all[breakpoint_all$largeall == "nolargeall",])
  samplesv_nolargeall = cbind(samplesv_nolargeall,in_out(samplesv_nolargeall))
 
  largeall_in <- nrow(samplesv_largeall[samplesv_largeall$chipseq >0,])
  largeall_out <- nrow(samplesv_largeall[samplesv_largeall$chipseq ==0,])
  nolargeall_in <- nrow(samplesv_nolargeall[samplesv_nolargeall$chipseq >0,])
  nolargeall_out <- nrow(samplesv_nolargeall[samplesv_nolargeall$chipseq ==0,])
  
  result <- fisher.test(matrix(data=c(largeall_in,
                                      largeall_out,
                                      nolargeall_in,
                                      nolargeall_out),,nrow = 2,ncol = 2))
                               
  if (result$p.value < 0.05 & result$estimate > 1) {
    resultadd <- data.frame("largeall_in"=largeall_in,
                            "largeall_out"=largeall_out,
                            "nolargeall_in"=nolargeall_in,
                            "nolargeall_out"=nolargeall_out,
                            "p.value"=result$p.value,
                            "or"=result$estimate,
                            "file"=filei)
    resultout <- rbind(resultout,resultadd)
  }
  
}




 ## ================== Load mappable file =================
# should first find out large SV breakpoint gathering point
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

checkNinSeq <- function(simple_sv){
  chr <- paste0("chr",simple_sv["chrom1"])
  pos1 <- min(simple_sv["pos1"],simple_sv["pos2"])
  pos2 <- max(simple_sv["pos1"],simple_sv["pos2"])
  cat(chr)
  out = unmasked(genome[[chr]])[pos1:pos2]
  print(out)
}

for (bioidi in bioid[1]) {
  simple_sv <- read.delim(file.path(root_dir,"scratch","sv-largeall",paste0(bioidi,"_largeall.tsv")))
  simple_sv <- simple_sv[simple_sv$SVtype %in% c("DEL","DUP","h2hINV","t2tINV"),]
  apply(simple_sv, 1, checkNinSeq)
  
}
