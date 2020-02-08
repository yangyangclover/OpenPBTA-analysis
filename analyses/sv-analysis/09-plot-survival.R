library(ggplot2)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(dplyr)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_directory <- file.path(root_dir, "scratch","sv-survival")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


## ===================== Load SV signatures =====================
sv_signature_path <- file.path("D:","Yang","Downloads","SigProfiler_2_5_1_7","SigProfiler_2_5_1_7","output","PBTA_PASS_v7.112","text","res_PBTA_PASS_v7.112_signature_activities_for_5_sigs.csv")
sv_signature <- read.csv(sv_signature_path)
sv_signature$Sample.Names <- paste0("BS_",sv_signature$Sample.Names)
row.names(sv_signature) <- sv_signature$Sample.Names
sv_signature <-t(sv_signature[,-1,drop=FALSE])


## ===================== Cluster:consensusClass =====================
title <- "sv signatures clusters"
results <- ConsensusClusterPlus(sv_signature,maxK = 7,reps = 50,pItem=0.8,pFeature = 1,
                                title = title, clusterAlg = "hc",distance = "pearson", seed = 123456,plot = "png")
clusternumber <- 5
sv_signature <- rbind(sv_signature,results[[clusternumber]][["consensusClass"]])
row.names(sv_signature)[nrow(sv_signature)] <- "consensusClass"


## ===================== Load Informations =====================
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
# creat info_histology
info_histology <- t(info[info$Kids_First_Biospecimen_ID %in% colnames(sv_signature),c("Kids_First_Biospecimen_ID","short_histology","age_at_diagnosis_days","germline_sex_estimate","OS_days","OS_status")])
colnames(info_histology) <- info_histology["Kids_First_Biospecimen_ID",]
# order info_histology
info_histology<-info_histology[,order(colnames(info_histology))]
# creat sv_signature_info
sv_signature_info  <- rbind(sv_signature,info_histology)
sv_signature_info <- t(sv_signature_info)
sv_signature_info <- as.data.frame(sv_signature_info)
# OS_status LIVING > 1, DIE >  2
sv_signature_info$OS_status <- gsub('LIVING', 1, sv_signature_info$OS_status)
sv_signature_info$OS_status <- gsub('DECEASED', 2, sv_signature_info$OS_status)
# read snv_df
snv_df <- read.csv(file.path(root_dir, "scratch","sv-snv","snv_in_samples_and_genes.csv"),row.names=1)
# sum all snv
sample_sumsv <- as.data.frame(colSums(snv_df))
# add  sample_sumsv
sv_signature_info <- merge(sv_signature_info,sample_sumsv,by = 0,all = TRUE)
sv_signature_info[is.na(sv_signature_info$`colSums(snv_df)`),"colSums(snv_df)"] <- 0
# add mutation_rate
sv_signature_info$mutation_rate <- sv_signature_info$`colSums(snv_df)`/30.88286401
# add short_histology_simple
sv_signature_info$short_histology_simple <- sv_signature_info$short_histology
sv_signature_info[as.character(sv_signature_info[,"short_histology_simple"]) %in% names(table(sv_signature_info[,"short_histology"])[order(table(sv_signature_info[,"short_histology"]),decreasing=T)[1:10]]) ==  FALSE,"short_histology_simple"] <- "Other"
# add disease simple
for (i in 1:length(unique(sv_signature_info$short_histology_simple))) {
  disease  <-  as.character(unique(sv_signature_info$short_histology_simple)[i])
  sv_signature_info[,disease] <- 0
  sv_signature_info[sv_signature_info[,"short_histology_simple"] == disease,disease]<-1
}
# add cluster
for (i in 1:length(unique(sv_signature_info$consensusClass))) {
  cluster  <-  paste0("cluster",as.character(unique(sv_signature_info$consensusClass)[i]))
  sv_signature_info[,cluster] <- 0
  sv_signature_info[sv_signature_info[,"consensusClass"] == as.character(unique(sv_signature_info$consensusClass)[i]),cluster]<-1
}
# add age_year
sv_signature_info$age_at_diagnosis_years <- floor(as.numeric(as.character(sv_signature_info$age_at_diagnosis_days))/365)
# change sex code
sv_signature_info$germline_sex_estimate <- sub("Female", 0, sv_signature_info$germline_sex_estimate) 
sv_signature_info$germline_sex_estimate <- sub("Male", 1, sv_signature_info$germline_sex_estimate)
# save sv_signature_info
write.csv(sv_signature_info,file.path(output_directory,"sv_signature_info.csv"),quote = F)


## ===================== Analysis =====================
Survdata <- Surv(as.numeric(as.character(sv_signature_info$OS_days)),as.numeric(sv_signature_info$OS_status))

sv_signature_info_disease <- sv_signature_info[sv_signature_info$short_histology != "LGAT" &sv_signature_info$short_histology != "HGAT",]
Survdata <- Surv(as.numeric(as.character(sv_signature_info_disease$OS_days)),as.numeric(sv_signature_info_disease$OS_status))
fit_KM <- survfit(Survdata~sv_signature_info_disease$cluster5,data=sv_signature_info_disease)

#[KM]
fit_KM <- survfit(Survdata~sv_signature_info$consensusClass,data=sv_signature_info)

#[Cox]
fit_Cox <- coxph(Survdata~
               #sv_signature_info$cluster1 +
               sv_signature_info$cluster5 +
               sv_signature_info$cluster4 +
               sv_signature_info$cluster3 +
               sv_signature_info$cluster2 +
               sv_signature_info$Ependymoma +
               sv_signature_info$LGAT +
               sv_signature_info$HGAT +
               sv_signature_info$Medulloblastoma +
               sv_signature_info$ATRT +
               sv_signature_info$Craniopharyngioma +
               sv_signature_info$Ganglioglioma +
               sv_signature_info$DNET +
               sv_signature_info$Meningioma +
               # sv_signature_info$Other +
               sv_signature_info$mutation_rate+
               sv_signature_info$age_at_diagnosis_years+
               sv_signature_info$germline_sex_estimate
             , 
             data=sv_signature_info)



#[Cox-Disease]
disease_simple <- unique(sv_signature_info[,"short_histology_simple"])
for  (diseasei  in  disease_simple) {
  DISEASE <- sv_signature_info[sv_signature_info[,diseasei] == 1,]
  DISEASE_Survdata <- Surv(as.numeric(as.character(DISEASE$OS_days)),as.numeric(DISEASE$OS_status))
  DISEASE_fit_KM <- survfit(DISEASE_Survdata~consensusClass,data=DISEASE)
  disease_survial_plot <- ggsurvplot(DISEASE_fit_KM, data = DISEASE,
             surv.median.line = "hv", # Add medians survival
             # Change legends: title & labels
             legend.title = "Signature Clusters",
             #legend.labs = c("cluster1", "cluster2","cluster3", "cluster4", "cluster5"),
             # Add p-value and tervals
             pval = TRUE,pval.size = 3,
             # Change censor
             censor.shape = 124,censor.size = 2,
             conf.int = FALSE,# 有无置信区间
             # break.x.by = 100, #横轴坐标
             # Add risk table
             risk.table = TRUE,tables.height = 0.3,tables.theme = theme_cleantable(),
             #palette = c("#E7B800", "#2E9FDF"),
             palette = c("red", "darkgreen","blue","deepskyblue","gold"),
             ggtheme = theme_bw(), # Change ggplot2 theme
             # Change font size, style and color
             main = "Survival curve",
             font.main = c(16, "bold", "darkblue"),
             font.x = c(14, "bold", "black"),
             font.y = c(14, "bold", "black"),
             font.tickslab = c(12, "plain", "black")
  )
  ggsave(file = paste0(diseasei,"_","PBTA_signature_5clusters_survival.png"), print(disease_survial_plot),path = output_directory,width = 15,height = 10)
  
  DISEASE_fit_KM_cluster5 <- survfit(DISEASE_Survdata~cluster5,data=DISEASE)
  disease_survial_plot <- ggsurvplot(DISEASE_fit_KM_cluster5, data = DISEASE,
                                     surv.median.line = "hv", # Add medians survival
                                     # Change legends: title & labels
                                     legend.title = "Signature Clusters",
                                     #legend.labs = c("cluster1", "cluster2","cluster3", "cluster4", "cluster5"),
                                     # Add p-value and tervals
                                     pval = TRUE,pval.size = 3,
                                     # Change censor
                                     censor.shape = 124,censor.size = 2,
                                     conf.int = FALSE,# 有无置信区间
                                     # break.x.by = 100, #横轴坐标
                                     # Add risk table
                                     risk.table = TRUE,tables.height = 0.3,tables.theme = theme_cleantable(),
                                     #palette = c("#E7B800", "#2E9FDF"),
                                     palette = c("red", "darkgreen","blue","deepskyblue","gold"),
                                     ggtheme = theme_bw(), # Change ggplot2 theme
                                     # Change font size, style and color
                                     main = "Survival curve",
                                     font.main = c(16, "bold", "darkblue"),
                                     font.x = c(14, "bold", "black"),
                                     font.y = c(14, "bold", "black"),
                                     font.tickslab = c(12, "plain", "black")
  )
  ggsave(file = paste0(diseasei,"_","PBTA_signature_clusters5_survival.png"), print(disease_survial_plot),path = output_directory,width = 15,height = 10)
  
  DISEASE_fit_Cox <- coxph(DISEASE_Survdata~
                             #sv_signature_info$cluster1 +
                             DISEASE$cluster5 +
                             DISEASE$cluster4 +
                             DISEASE$cluster3 +
                             DISEASE$cluster2 +
                             DISEASE$mutation_rate+
                             DISEASE$age_at_diagnosis_years+
                             DISEASE$germline_sex_estimate, 
                           data=DISEASE)
  
  # Prepare the columns
  HR <- round(exp(coef(DISEASE_fit_Cox)), 2)
  CI <- round(exp(confint(DISEASE_fit_Cox)), 2)
  P <- round(coef(summary(DISEASE_fit_Cox))[,5], 3)
  # Names the columns of CI
  colnames(CI) <- c("Lower", "Higher")
  # Bind columns together as dataset
  table2 <- as.data.frame(cbind(HR, CI, P))
  write.table(table2,file=file.path(output_directory,paste0(diseasei,"_fit_Cox.txt")))
}


## ===================== Plot Survial Curve =====================
survial_plot  <- ggsurvplot(fit_KM, data = sv_signature_info_disease,
           surv.median.line = "hv", # Add medians survival
           # Change legends: title & labels
           legend.title = "Signature Clusters",
           #legend.labs = c("cluster1", "cluster2","cluster3", "cluster4", "cluster5"),
           # Add p-value and tervals
           pval = TRUE,pval.size = 3,
           # Change censor
           censor.shape = 124,censor.size = 2,
           conf.int = FALSE,# 有无置信区间
           # break.x.by = 4, #横轴坐便
           # Add risk table
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           #palette = c("#E7B800", "#2E9FDF"),
           palette = c("firebrick3","coral","dodgerblue3","darkorchid","darkgreen"),
           ggtheme = theme_bw(), # Change ggplot2 theme
           # Change font size, style and color
           main = "Survival curve",
           font.main = c(16, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black")
)
ggsave(file = "PBTA_signature_5clusters_survival.jpeg", print(survial_plot),path = output_directory,width = 12,height = 10)





XX<-cor(sv_signature_info[,c("cluster5","HGAT")])
XX<-cor(sv_signature_info[,c("age_at_diagnosis_years","mutation_rate")])
XX<-cor(sv_signature_info[,c(14,16:31)])
kappa(XX,exact=TRUE) 
heatmap(XX)
