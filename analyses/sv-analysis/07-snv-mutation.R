library(rprojroot)
library(readr)
library(magrittr)
library(dplyr)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(ggplot2)
library(randomcoloR)
library(RColorBrewer)
library(circlize)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_directory <- file.path(root_dir, "scratch","sv-snv")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data","independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load SNV File and Generate snv_df =====================
# You can skip this section if you have already had a snv_file in output_directory
# OK tell me what kind of SNV mutation you like
# you can choose all or driver
# write it down
snvtype <-  "all" # you can write down driver or all

# load snv maf file
snv_file <- read_tsv(file.path(root_dir,"data", "pbta-snv-consensus-mutation.maf.tsv.gz")) %>% filter(Tumor_Sample_Barcode  %in% bioid)

# if you write down "driver", if you don't, run it also ok
if (snvtype == "driver") {
  driver_mutation_type <- c("Frame_Shift_Del","Frame_Shift_Ins",
                            "Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region",
                            "Splice_Site","Translation_Start_Site")
  snv_file_driver <- snv_file[snv_file$Variant_Classification %in% driver_mutation_type,]
  snv_file <- snv_file_driver
  csvname  <- "driver_snv_in_samples_and_genes.csv"
} else {
  csvname <- "snv_in_samples_and_genes.csv"
}

# creat snv_df, row is gene, col is sample
snv_df <- data.frame(matrix(0,ncol = length(unique(snv_file$Tumor_Sample_Barcode)),nrow=length(unique(snv_file$Hugo_Symbol))))
row.names(snv_df) <- unique(snv_file$Hugo_Symbol)
names(snv_df) <-unique(snv_file$Tumor_Sample_Barcode)

# count gene snv for all samples
for (i in 1:nrow(snv_file)) {
  sample  <- as.character(snv_file[i,"Tumor_Sample_Barcode"])
  gene <- as.character(snv_file[i,"Hugo_Symbol"])
  snv_df[gene,sample] <-  snv_df[gene,sample] +  1
}

# write csv file
write.csv(snv_df,file.path(output_directory,csvname),quote=F)


## ================== Make a new snv_df_new, with sum information ====================
# if you have already had a snv_df csv file in your output_directory, let's read it
# OK tell me what kind of SNV mutation you like
# you can choose all or driver
# write it down
snvtype <-  "driver" # you can write down driver or all
if (snvtype == "driver") {
  snv_df<-read.csv(file.path(output_directory,"driver_snv_in_samples_and_genes.csv"),row.names=1)
} else {
  snv_df<-read.csv(file.path(output_directory,"snv_in_samples_and_genes.csv"),row.names=1)
}
# creat snv_df_new
snv_df_new <- snv_df
# sum col and row
snv_df_new$sum <- rowSums(snv_df_new[,c(1:ncol(snv_df_new))])
snv_df_new["sum",c(1:ncol(snv_df_new))] <-  colSums(snv_df_new[,c(1:ncol(snv_df_new))])
# sum col and row only if value > 0
snv_df_new$uniquesum <- rowSums(snv_df_new[,c(1:ncol(snv_df_new)-1)] > 0)
snv_df_new["uniquesum",c(1:ncol(snv_df_new))] <-colSums(snv_df_new[c(1:nrow(snv_df_new)-1),c(1:ncol(snv_df_new))] > 0)


## ===================== Load SV signatures =====================
sv_signature_path <- file.path("D:","Yang","Downloads","SigProfiler_2_5_1_7","SigProfiler_2_5_1_7","output","PBTA_PASS_v7.112","text","res_PBTA_PASS_v7.112_signature_activities_for_5_sigs.csv")
sv_signature <- read.csv(sv_signature_path)
sv_signature$Sample.Names <- paste0("BS_",sv_signature$Sample.Names)
row.names(sv_signature) <- sv_signature$Sample.Names
sv_signature <-t(sv_signature[,-1,drop=FALSE])


## ===================== Cluster:consensusClass =====================
title <- "SV Signatures Clusters"
# do cluster
results <- ConsensusClusterPlus(sv_signature,maxK = 7,reps = 50,pItem=0.8,pFeature = 1,
                                title = title, clusterAlg = "hc",distance = "pearson", seed = 123456,plot = "png")

#  choose cluster number
clusternumber <- 5

#  add cluster result to sv_signature
sv_signature <- rbind(sv_signature,results[[clusternumber]][["consensusClass"]])
row.names(sv_signature)[nrow(sv_signature)] <- "consensusClass"

# # add  name
# sv_signature <- rbind(colnames(sv_signature),sv_signature)
# rownames(sv_signature)[1] <- "Name"

## ===================== Combine SNV and SV， get sv_snv =====================
# snv_df_simple: value > 1 ===>>>>>  value = 1
snv_df_simple <- snv_df
snv_df_simple[snv_df_simple>=1]<-1
# creat sv_snv
sv_snv <- bind_rows(as.data.frame(sv_signature), snv_df_simple,snv_df_new["sum",c(1:(ncol(snv_df_new)-2))])
row.names(sv_snv) <- c(row.names(sv_signature),row.names(snv_df_simple),"sum")
# since bind_row will generate NA, so NA ====>>>> 0
sv_snv[is.na(sv_snv)] <- 0
# add mutation_rate
sv_snv["mutation_rate",] <- sv_snv["sum",]/30.88286401
# add mutation_rate_log10
sv_snv["mutation_rate_log10",] <- log10(as.numeric(sv_snv["mutation_rate",])+1)


## ===================== Do Fisher Test =====================
# If you have had a Fisher Test result in your output_directory, skip this section
# output file in this section will not be  used in the following sections
# Do Fisher Test
potentical_gene <- row.names(snv_df)
# significant_association is used to association resutl
significant_association <- data.frame()
for (classi in 1:clusternumber)  {
  loop <- 0
  for (genei in potentical_gene ){
    loop <- loop +1
    gene1_cluster1  <- sum((sv_snv["consensusClass",] == classi) & (sv_snv[genei,] == 1))
    gene0_cluster1  <- sum((sv_snv["consensusClass",] == classi) & (sv_snv[genei,] == 0))
    gene1_cluster0  <- sum((sv_snv["consensusClass",] != classi) & (sv_snv[genei,] == 1))
    gene0_cluster0  <- sum((sv_snv["consensusClass",] != classi) & (sv_snv[genei,] == 0))
    testmatrix <- c(gene1_cluster1,gene0_cluster1,gene1_cluster0,gene0_cluster0)
    dim(testmatrix) <- c(2,2)
    test <-fisher.test(testmatrix)
    p <-test[1]
    or <- test[3]
    significant_association <- rbind(
      significant_association,
      data.frame(
        "gene" = genei,
        "cluster" = "consensusClass",
        "class" = classi,
        "p" = p,
        "OR" = or
      )
    )
  }
}
# log10(p)
significant_association[,"log10p"] <- -log10(significant_association$p.value)
# Fdr(p) in individual cluster
significant_association[significant_association$class == 1,"FDR"] <- p.adjust(significant_association[significant_association$class == 1,"p.value"],method = "fdr",n=nrow(significant_association[significant_association$class == 1,]))
significant_association[significant_association$class == 2,"FDR"] <- p.adjust(significant_association[significant_association$class == 2,"p.value"],method = "fdr",n=nrow(significant_association[significant_association$class == 2,]))
significant_association[significant_association$class == 3,"FDR"] <- p.adjust(significant_association[significant_association$class == 3,"p.value"],method = "fdr",n=nrow(significant_association[significant_association$class == 3,]))
significant_association[significant_association$class == 4,"FDR"] <- p.adjust(significant_association[significant_association$class == 4,"p.value"],method = "fdr",n=nrow(significant_association[significant_association$class == 4,]))
significant_association[significant_association$class == 5,"FDR"] <- p.adjust(significant_association[significant_association$class == 5,"p.value"],method = "fdr",n=nrow(significant_association[significant_association$class == 5,]))
# write csv
snvtype <-  "driver" # you can write down driver or all
if (snvtype == "driver") {
  csvname <- "significant_association_value_fishertest_drivermutation.csv"
} else {
  csvname <- "significant_association_value_fishertest.csv"
}
write.csv(significant_association,file.path(output_directory,csvname),quote=F)


## ===================== Load Informations =====================
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
info_histology <- t(info[info$Kids_First_Biospecimen_ID %in% colnames(sv_snv),c("Kids_First_Biospecimen_ID","short_histology","age_at_diagnosis_days")])
# add age_at_diagnosis_years
info_histology <- rbind(info_histology,as.numeric(info_histology["age_at_diagnosis_days",])/365)
rownames(info_histology)[4] <- "age_at_diagnosis_years"
#info_histology["age_at_diagnosis_years",as.numeric(info_histology["age_at_diagnosis_years",])>18] <- 18
# add short_histology_simple
# top 10 will be left, others will be other
info_histology <- rbind(info_histology,info_histology["short_histology",])
rownames(info_histology)[5] <- "short_histology_simple"
info_histology["short_histology_simple",as.character(info_histology["short_histology_simple",]) %in% names(table(info_histology["short_histology",])[order(table(info_histology["short_histology",]),decreasing=T)[1:10]]) ==  FALSE] <- "Other"
# set colname
colnames(info_histology) <- info_histology[1,]
# order sv_signature
sv_signature <- sv_signature[,order(colnames(sv_signature))]
# order info_histology
info_histology<-info_histology[,order(colnames(info_histology))]

## ===================== DO test for Disease v.s Cluster =====================
# If you have had a Fisher Test result in your output_directory, skip this section
# output file in this section will not be  used in the following sections
# Do Fisher Test
# creat disease_cluster
disease_cluster <- rbind(sv_signature,info_histology)
rownames(disease_cluster) <-  c("consensusClass","short_histology")
# potentical_disease is all tumor types
potentical_disease <- unique(disease_cluster["short_histology",])
# disease_cluster_association is used to save association result
disease_cluster_association <- data.frame()
for (classi in 1:clusternumber)  {
  loop <- 0
  for (diseasei in potentical_disease ){
    loop <- loop +1
    disease1_cluster1  <- sum((disease_cluster["consensusClass",] == classi) & (disease_cluster["short_histology",] == as.character(diseasei)))
    disease0_cluster1  <- sum((disease_cluster["consensusClass",] == classi) & (disease_cluster["short_histology",] != as.character(diseasei)))
    disease1_cluster0  <- sum((disease_cluster["consensusClass",] != classi) & (disease_cluster["short_histology",] == as.character(diseasei)))
    disease0_cluster0  <- sum((disease_cluster["consensusClass",] != classi) & (disease_cluster["short_histology",] != as.character(diseasei)))
    testmatrix <- c(disease1_cluster1,disease0_cluster1,disease1_cluster0,disease0_cluster0)
    dim(testmatrix) <- c(2,2)
    test <-fisher.test(testmatrix)
    p <-test[1]
    or <- test[3]
    disease_cluster_association <- rbind(
      disease_cluster_association,
      data.frame(
        "disease" = diseasei,
        "cluster" = "consensusClass",
        "class" = classi,
        "disease1_cluster1" =disease1_cluster1,
        "disease0_cluster1" =disease0_cluster1,
        "disease1_cluster0" =disease1_cluster0,
        "disease0_cluster0" =disease0_cluster0,
        "p" = p,
        "OR" = or
      )
    )
  }
}
# log10(p)
disease_cluster_association[,"log10p"] <- -log10(disease_cluster_association$p.value)
# fdr(p)
disease_cluster_association[disease_cluster_association$class == 1,"FDR"] <- p.adjust(disease_cluster_association[disease_cluster_association$class == 1,"p.value"],method = "fdr",n=nrow(disease_cluster_association[disease_cluster_association$class == 1,]))
disease_cluster_association[disease_cluster_association$class == 2,"FDR"] <- p.adjust(disease_cluster_association[disease_cluster_association$class == 2,"p.value"],method = "fdr",n=nrow(disease_cluster_association[disease_cluster_association$class == 2,]))
disease_cluster_association[disease_cluster_association$class == 3,"FDR"] <- p.adjust(disease_cluster_association[disease_cluster_association$class == 3,"p.value"],method = "fdr",n=nrow(disease_cluster_association[disease_cluster_association$class == 3,]))
disease_cluster_association[disease_cluster_association$class == 4,"FDR"] <- p.adjust(disease_cluster_association[disease_cluster_association$class == 4,"p.value"],method = "fdr",n=nrow(disease_cluster_association[disease_cluster_association$class == 4,]))
disease_cluster_association[disease_cluster_association$class == 5,"FDR"] <- p.adjust(disease_cluster_association[disease_cluster_association$class == 5,"p.value"],method = "fdr",n=nrow(disease_cluster_association[disease_cluster_association$class == 5,]))
# write csv
write.csv(disease_cluster_association,file.path(output_directory,"disease_cluster_association_fishertest.csv"),quote=F)



## ===================== Plot heatmap ===================== 
# ------------significantgene------------
# Cluster4      AL590867.1
# Cluster4      AC068631.2
# Cluster4      TPI1P1
# Cluster4      HNRNPCP2
# Cluster4      CTSB
# Cluster4      RPS23P8
# Cluster4      BGN
# Cluster4      AC242426.2
# Cluster4      AL450998.1
# Cluster5      PIK3C2G
# Cluster5      LINC01781
# Cluster5      CNTNAP5
# Cluster5      PCSK5
# Cluster5      GRIN2A
# Cluster5      MIR99AHG
# Cluster5      TP53
# --------------------------------------
# Cluster5      H3F3A
# Cluster5      TP53
# -------------------------------------
sv_snv <- as.matrix(sv_snv)
sv_snv_cluster_info <- rbind(info_histology,sv_signature,sv_snv)


disease_col <- c("blue","darkorange","cornsilk","brown3","skyblue","limegreen","cadetblue","gold","deeppink","yellow")
column_ha <- HeatmapAnnotation("consensusClass"=sv_snv["consensusClass",],
                               # "AL590867.1"= sv_snv["AL590867.1",],
                               # "AC068631.2"= sv_snv["AC068631.2",],
                               # "TPI1P1"= sv_snv["TPI1P1",],
                               # "HNRNPCP2"= sv_snv["HNRNPCP2",],
                               # "CTSB"= sv_snv["CTSB",],
                               # "RPS23P8"= sv_snv["RPS23P8",],
                               # "BGN"= sv_snv["BGN",],
                               # "AC084880.2"= sv_snv["AC084880.2",],
                               # "AL450998.1"= sv_snv["AL450998.1",],
                               # "PIK3C2G"= sv_snv["PIK3C2G",],
                               # "LINC01781"= sv_snv["LINC01781",],
                               # "CNTNAP5"= sv_snv["CNTNAP5",],
                               # "PCSK5"= sv_snv["PCSK5",],
                               # "GRIN2A"= sv_snv["GRIN2A",],
                               # "MIR99AHG"= sv_snv["MIR99AHG",],
                               "BRCA1"= sv_snv["BRCA1",],
                               "BRCA2"= sv_snv["BRCA2",],
                               "ATRX"= sv_snv["ATRX",],
                               "H3F3A"= sv_snv["H3F3A",],
                               "PARP1"= sv_snv["TP53",],
                               "TP53"= sv_snv["PARP1",],
                               # "BRCA1"= sv_snv["BRCA1",],
                               # "BRCA2"= sv_snv["BRCA2",],
                               # "TTN"= sv_snv["TTN",],
                               # "PCLO"= sv_snv["PCLO",],
                               # "IGFBP7"= sv_snv["IGFBP7",],
                               "disease" = info_histology["short_histology_simple",],
                               "age" = as.numeric(info_histology["age_at_diagnosis_years",]),
                               "log10(mutation_rate)"  = anno_barplot(log10(as.numeric(sv_snv["sum",])+0.1)),
                               col = list(consensusClass = c("1" = "firebrick3", "2" = "coral", "3" = "dodgerblue3",
                                                             "4" = "darkorchid", "5" = "darkgreen"),
                                          # AL590867.1 = c("0"="white","1"="black"),
                                          # AC068631.2 = c("0"="white","1"="black"),
                                          # TPI1P1 = c("0"="white","1"="black"),
                                          # HNRNPCP2= c("0"="white","1"="black"),
                                          # CTSB = c("0"="white","1"="black"),
                                          # RPS23P8= c("0"="white","1"="black"),
                                          # BGN= c("0"="white","1"="black"),
                                          # AC084880.2= c("0"="white","1"="black"),
                                          # AL450998.1= c("0"="white","1"="black"),
                                          # PIK3C2G = c("0"="white","1"="black"),
                                          # LINC01781= c("0"="white","1"="black"),
                                          # CNTNAP5= c("0"="white","1"="black"),
                                          # PCSK5= c("0"="white","1"="black"),
                                          # GRIN2A= c("0"="white","1"="black"),
                                          # MIR99AHG= c("0"="white","1"="black"),
                                          TP53= c("0"="white","1"="black"),
                                          H3F3A= c("0"="white","1"="black"),
                                          BRCA1= c("0"="white","1"="black"),
                                          BRCA2= c("0"="white","1"="black"),
                                          PARP1= c("0"="white","1"="black"),
                                         # TTN= c("0"="white","1"="black"),
                                         # PCLO= c("0"="white","1"="black"),
                                         # IGFBP7= c("0"="white","1"="black"),
                                         ATRX= c("0"="white","1"="black"),
                                          disease = setNames(disease_col, unique(info_histology["short_histology_simple",]))
                                          ),
                               show_legend = c(rep(FALSE,7),rep(TRUE,4))#,
                               # annotation_legend_param = list(disease = list(
                               #   title = "disease",
                               #   at = disease_col,
                               #   labels = unique(info_histology["short_histology_simple",])
                               # ))
)

column_his <- HeatmapAnnotation(
  sig_number = anno_barplot(t(sv_snv[c(1:5),]),
                            height = unit(5, "cm"),
                            width = 1,
                            gp = gpar(col=c("firebrick3","dodgerblue3","coral","darkorchid","darkgreen"))),
  show_legend = TRUE
)

zero_row_mat <- matrix(nrow = 0, ncol = ncol(sv_snv))

ht <- Heatmap(zero_row_mat,
              top_annotation = column_ha,
              bottom_annotation = column_his,
        column_order = results[[clusternumber]]$consensusTree$order,
        show_heatmap_legend = FALSE,
        show_column_names = FALSE,
        clustering_distance_columns = "pearson" # ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
        )
draw(ht,
     show_annotation_legend = TRUE,
     annotation_legend_side = "bottom")

