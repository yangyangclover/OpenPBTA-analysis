setwd("~/GitHub/OpenPBTA-analysis")

library(ConsensusClusterPlus)

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


## ===================== Read simple SV and Keep Large Del/Dup/Inv===================
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

## ==================== Load Signature and clusters ====================
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


## =================== Load H3 mutation ==================
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

## ================= Load Disease Info ==================
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
# combine disease
largeall_stat <- merge(largeall_stat,info[,c("Kids_First_Biospecimen_ID","short_histology")],by.x = "bioid",by.y = "Kids_First_Biospecimen_ID",all.x = T)


write.csv(largeall_stat,file.path(output_directory,"patients_largeall_cluster_H3mutation_disease.csv"),quote = F,row.names = F)
