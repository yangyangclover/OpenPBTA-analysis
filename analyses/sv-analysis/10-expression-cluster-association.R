library(plyr)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# load file: patientid2wgsid
sample_list <- read_tsv(file.path(root_dir,"data","independent-specimens.wgs.primary-plus.tsv"))
# load file: allid2patientid
sample_info <-  read_tsv(file.path(root_dir,"data","pbta-histologies.tsv"))
# load file: gene expression data
expression <- readRDS(file.path(root_dir,"data","pbta-gene-expression-rsem-fpkm.stranded.rds"))

# creat idlab
rna_id <- colnames(expression)[2:ncol(expression)]
idlab <- as.data.frame(rna_id)
idlab$patient_id <- sample_info[sample_info$Kids_First_Biospecimen_ID %in% rna_id,]$Kids_First_Participant_ID
idlab <- merge(idlab,sample_list,all = FALSE, by.x="patient_id", by.y="Kids_First_Participant_ID")
idlab <- merge(idlab,sample_info[,c("Kids_First_Biospecimen_ID","tumor_descriptor")],by.x="rna_id",by.y="Kids_First_Biospecimen_ID",all = FALSE,)
colnames(idlab)[4]<-"rna_tumor_descriptor"
idlab <- merge(idlab,sample_info[,c("Kids_First_Biospecimen_ID","tumor_descriptor")],by.x="Kids_First_Biospecimen_ID",by.y="Kids_First_Biospecimen_ID",all = FALSE,)
colnames(idlab)[5]<-"wgs_tumor_descriptor"

# duplicate wgs
rep <- idlab[idlab$Kids_First_Biospecimen_ID %in% count(idlab$Kids_First_Biospecimen_ID)[count(idlab$Kids_First_Biospecimen_ID)$freq>1,"x"],]
# idlab_unique
idlab_unique  <- data.frame()
for (i in unique(idlab$Kids_First_Biospecimen_ID)) {
 if (i %in% rep$Kids_First_Biospecimen_ID == FALSE) {
   add <- idlab[idlab$Kids_First_Biospecimen_ID == i,]
   idlab_unique <- rbind(idlab_unique,add)
 } else {
   if  (sum(rep[rep$Kids_First_Biospecimen_ID == i,"rna_tumor_descriptor"] == rep[rep$Kids_First_Biospecimen_ID == i, "wgs_tumor_descriptor"]) > 0 ){
     add <-  rep[rep$Kids_First_Biospecimen_ID == i,]
     add <- add[add$rna_tumor_descriptor == add$wgs_tumor_descriptor,]
     idlab_unique <- rbind(idlab_unique,add[sample(1:nrow(add),size = 1),])
   }  else {
     add  <- rep[rep$Kids_First_Biospecimen_ID == i,]
     idlab_unique <- rbind(idlab_unique,add[sample(1:nrow(add),size = 1),])
   }
 }
}

# load cluster
cluster <- read.csv(file = file.path(root_dir,"scratch","sv-survival","sv_signature_info.csv"))

#merge idlab_unique and cluster
idlab_unique_cluster <- merge(idlab_unique,cluster[,c("Kids_First_Biospecimen_ID","consensusClass","Signature.1","Signature.2","Signature.3","Signature.4","Signature.5","short_histology")],by.x="Kids_First_Biospecimen_ID",by.y="Kids_First_Biospecimen_ID")
idlab_unique_cluster <- t(idlab_unique_cluster)
colnames(idlab_unique_cluster) <- idlab_unique_cluster["Kids_First_Biospecimen_ID",]

#merge idlab_unique_cluster and expression
gene_expression <- expression[,c(2:ncol(expression))]
rownames(gene_expression) <- expression[,"gene_id"]

# Input your Gene
gene <-"PGBD5"

analysis <-  t(gene_expression[grep(pattern = gene,rownames(gene_expression),),])
analysis <- cbind(analysis,rownames(analysis))
colnames(analysis)[2] <-  "rna_id"
rnaid_cluster <- t(as.data.frame(idlab_unique_cluster)[c("rna_id","consensusClass"),])
rnaid_tumortype <- t(as.data.frame(idlab_unique_cluster)[c("rna_id","short_histology"),])

analysis_table<-merge(analysis,rnaid_cluster)
colnames(analysis_table)[2] <- "mygene"
analysis_table$mygene <- as.numeric(as.character(analysis_table[,2]))
analysis_table$cluster1 <- as.character(analysis_table$consensusClass)
analysis_table$cluster2 <- as.character(analysis_table$consensusClass)
analysis_table$cluster3 <- as.character(analysis_table$consensusClass)
analysis_table$cluster4 <- as.character(analysis_table$consensusClass)
analysis_table$cluster5 <- as.character(analysis_table$consensusClass)
analysis_table$cluster1[which(analysis_table$cluster2 !=1)] <- "uncluster1"
analysis_table$cluster1[which(analysis_table$cluster2 ==1)] <- "cluster1"
analysis_table$cluster2[which(analysis_table$cluster2 !=2)] <- "uncluster2"
analysis_table$cluster2[which(analysis_table$cluster2 ==2)] <- "cluster2"
analysis_table$cluster3[which(analysis_table$cluster2 !=3)] <- "uncluster3"
analysis_table$cluster3[which(analysis_table$cluster2 ==3)] <- "cluster3"
analysis_table$cluster4[which(analysis_table$cluster2 !=4)] <- "uncluster4"
analysis_table$cluster4[which(analysis_table$cluster2 ==4)] <- "cluster4"
analysis_table$cluster5[which(analysis_table$cluster2 !=5)] <- "uncluster2"
analysis_table$cluster5[which(analysis_table$cluster2 ==5)] <- "cluster2"

# Input my tumor type
mytumortype <- c("ATRT","Medulloblastoma")
analysis_table_tumortype<-merge(analysis,rnaid_tumortype)
colnames(analysis_table_tumortype)[2] <- "mygene"
analysis_table_tumortype$mygene <- as.numeric(as.character(analysis_table_tumortype[,2]))
analysis_table_tumortype$mytumortype <- as.character(analysis_table_tumortype$short_histology)
analysis_table_tumortype$mytumortype[which(analysis_table_tumortype$mytumortype %in% mytumortype == FALSE )] <- "Others"


ggplot(data=analysis_table_tumortype,aes(x=mytumortype, y=mygene)) + geom_violin(aes(fill=factor(mytumortype)))+ 
  geom_violin(aes(fill=factor(mytumortype))) + 
  labs(x = "Tumor Type", y = "Gene Expression", title=gene) +
  theme(axis.text.x  = element_text(angle=30, vjust=0.5) )

ggplot(data=analysis_table_tumortype,aes(x=mytumortype, y=mygene)) + 
  geom_boxplot(aes(fill=factor(mytumortype)))+ 
  labs(x = "Tumor Type", y = "Gene Expression", title=gene) +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5,size = 10),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        axis.line.x=element_line(colour="black",size = 1),
        axis.line.y=element_line(colour="black",size = 1),
        panel.background=element_blank(),
        legend.background = element_blank())+
  scale_fill_manual(values = list("ATRT"="firebrick3",
                                  "Medulloblastoma"="dodgerblue3",
                                  "Others"="darkgreen"))

ggplot(data=analysis_table_tumortype,aes(x=short_histology, y=mygene)) + 
  geom_violin(aes(fill=factor(short_histology))) + 
  labs(x = "Tumor Type", y = "Gene Expression", title=gene) +
  theme(axis.text.x  = element_text(angle=30, vjust=0.5)
        axis.line = element_text())

ggplot(data=analysis_table_tumortype,aes(x=short_histology, y=mygene)) + 
  geom_boxplot(aes(fill=factor(short_histology))) + 
  labs(x = "Tumor Type", y = "Gene Expression", title=gene) +
  theme(axis.text.x  = element_text(angle=30, vjust=0.5) )

TukeyHSD(aov(formula = mygene ~ short_histology, data = analysis_table_tumortype))
TukeyHSD(aov(formula = mygene ~ short_histology, data = analysis_table_tumortype))

# plot all clusters
ggplot(data=analysis_table,aes(x=consensusClass, y=mygene)) + geom_violin(aes(fill=factor(consensusClass)))
ggplot(data=analysis_table,aes(x=consensusClass, y=mygene)) + geom_boxplot(aes(fill=factor(consensusClass)))

# plot my cluster
ggplot(data=analysis_table,aes(x=cluster2, y=mygene)) + geom_violin(aes(fill=factor(cluster2)))
ggplot(data=analysis_table,aes(x=cluster2, y=mygene)) + 
  geom_boxplot(aes(fill=factor(cluster2))) +
  labs(x = "Cluster", y = "Gene Expression", title=gene) +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5,size = 20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        axis.line.x=element_line(colour="black",size=1),
        axis.line.y=element_line(colour="black",size=1),
        panel.background=element_blank(),
        legend.background = element_blank())+
  scale_fill_manual(values = list("cluster2"="firebrick3","uncluster2"="dodgerblue3"))


# ttest
# Input your cluster
testcluster <- 4
t.test(as.numeric(as.character(analysis_table[analysis_table$consensusClass ==testcluster,2])),
       as.numeric(as.character(analysis_table[analysis_table$consensusClass !=testcluster,2])))
TukeyHSD(aov(formula = mygene ~ consensusClass, data = analysis_table))

# =========================================My gene and sig number======================
gene <-"PGBD5"

rnaid_sig <- t(as.data.frame(idlab_unique_cluster)[c("rna_id","Signature.1","Signature.2","Signature.3","Signature.4","Signature.5"),])

gene_sig_analysis_table<-merge(analysis,rnaid_sig)
colnames(gene_sig_analysis_table)[2] <- "mygene"
gene_sig_analysis_table$mygene <- as.numeric(as.character(gene_sig_analysis_table$mygene))
gene_sig_analysis_table$Signature.1 <- as.numeric(as.character(gene_sig_analysis_table$Signature.1))
gene_sig_analysis_table$Signature.2 <- as.numeric(as.character(gene_sig_analysis_table$Signature.2))
gene_sig_analysis_table$Signature.3 <- as.numeric(as.character(gene_sig_analysis_table$Signature.3))
gene_sig_analysis_table$Signature.4 <- as.numeric(as.character(gene_sig_analysis_table$Signature.4))
gene_sig_analysis_table$Signature.5 <- as.numeric(as.character(gene_sig_analysis_table$Signature.5))
cor.test(gene_sig_analysis_table$mygene,gene_sig_analysis_table$Signature.3)
t.test(gene_sig_analysis_table[gene_sig_analysis_table$mygene>=0 & gene_sig_analysis_table$mygene<5.735,"Signature.3"],
         gene_sig_analysis_table[gene_sig_analysis_table$mygene>=5.735 & gene_sig_analysis_table$mygene<60,"Signature.3"])

# ============================================= age vs. cluster ======================
cluster_age  <- merge(cluster[,c("Kids_First_Biospecimen_ID","consensusClass","Signature.1","Signature.2","Signature.3","Signature.4","Signature.5","short_histology")],sample_info[,c("Kids_First_Biospecimen_ID","age_at_diagnosis_days")],by.x="Kids_First_Biospecimen_ID",by.y="Kids_First_Biospecimen_ID")
cluster_age$age_at_diagnosis_days <- as.numeric(cluster_age$age_at_diagnosis_days)
cluster_age$sigsum <- apply(cluster_age[,2:6],1,sum)

cor.test(cluster_age$sigsum,cluster_age$age_at_diagnosis_days, method=c("pearson", "kendall", "spearman"))
ggplot(data=cluster_age,aes(x=sigsum,y=age_at_diagnosis_days)) + geom_point()
ggplot(data=cluster_age,aes(x=consensusClass,y=age_at_diagnosis_days)) + geom_violin(aes(fill=factor(consensusClass)))
ggplot(data=cluster_age,aes(x=consensusClass,y=age_at_diagnosis_days)) + 
  geom_boxplot(aes(fill=factor(consensusClass)))  +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5,size = 20),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        axis.line.x=element_line(colour="black",size = 1),
        axis.line.y=element_line(colour="black",size = 1),
        panel.background=element_blank(),
        legend.background = element_blank())+
  scale_fill_manual(values = list("1"="firebrick3",
                                  "2"="coral",
                                  "3"="dodgerblue3",
                                  "4"="darkorchid",
                                  "5"="darkgreen"))

cor.test(cluster_age[cluster_age$consensusClass == 5,"Signature.5"], cluster_age[cluster_age$consensusClass == 5,"age_at_diagnosis_days"], method=c("pearson", "kendall", "spearman"))


cluster_age_melt <- melt(cluster_age, id=c("Kids_First_Biospecimen_ID","consensusClass","short_histology","age_at_diagnosis_days"))
colnames(cluster_age_melt)[5:6] <- c("signature","signature_number")
ggplot(cluster_age_melt, aes(x=signature_number, y=age_at_diagnosis_days)) + 
  geom_point() + 
  stat_smooth(method = "lm", se = TRUE) +
  facet_grid(signature ~ .) +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5,size = 15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.line.x=element_line(colour="black",size = 1),
        axis.line.y=element_line(colour="black",size = 1),
        panel.background=element_blank(),
        legend.background = element_blank())

# ============================= SV signature and mutation signature ===========================
mutsig <- read_tsv(file.path(root_dir,"analyses","mutational-signatures","results","cosmic_signatures_results.tsv"))
mutsig$signature <- paste0("Mutation_",mutsig$signature)
svsig <-  cluster_age_melt[,c("Kids_First_Biospecimen_ID","signature","signature_number")]
colnames(svsig) <- c("Tumor_Sample_Barcode","signature","num_mutations")
mutsig_svsig<-rbind(mutsig[,c("Tumor_Sample_Barcode","signature","num_mutations")],svsig)
mutsig_svsig_dcast <- dcast(mutsig_svsig,Tumor_Sample_Barcode ~  signature)
mutsig_svsig_dcast <- na.omit(mutsig_svsig_dcast)


mutsigidi <- 0
mutsig_svsig_ptable<-data.frame(matrix(NA,30,5))
mutsig_svsig_cortable<-data.frame(matrix(NA,30,5))
colnames(mutsig_svsig_ptable) <- colnames(mutsig_svsig_dcast[,32:36])
colnames(mutsig_svsig_cortable) <- colnames(mutsig_svsig_dcast[,32:36])
row.names(mutsig_svsig_ptable) <- colnames(mutsig_svsig_dcast[,2:31])
row.names(mutsig_svsig_cortable) <- colnames(mutsig_svsig_dcast[,2:31])
for (mutsigid  in mutsig_svsig_dcast[,2:31]) {
  mutsigidi  <- mutsigidi + 1
  svsigidi <-  0
  for (svsigid in mutsig_svsig_dcast[32:36]) {
    svsigidi <- svsigidi + 1
    cortest <- cor.test(mutsigid,svsigid,method=c("pearson"))
    cortestp <- cortest$p.value
    cortestcor <- cortest$estimate
    mutsig_svsig_ptable[mutsigidi,svsigidi] <- cortestp
    mutsig_svsig_cortable[mutsigidi,svsigidi] <- cortestcor
    if (cortestp < 0.05 & (cortestcor > 0.2|cortestcor < -0.2))  {
      print(colnames(mutsig_svsig_dcast[mutsigidi+1]))
      print(colnames(mutsig_svsig_dcast[svsigidi+31]))
      print(cortest)
    }
  }
}
write.csv(mutsig_svsig_ptable,"mutsig_svsig_ptable.csv",quote = FALSE)
write.csv(mutsig_svsig_cortable,"mutsig_svsig_cortable.csv",quote = FALSE)
ggplot(data=mutsig_svsig_dcast,aes(x=Mutation_Signature.13,y=Signature.1)) + geom_point()

mutsig_svsig_plot <- cbind(melt(mutsig_svsig_cortable ),melt(mutsig_svsig_ptable ))
mutsig_svsig_plot[,1] <- rep(row.names(mutsig_svsig_cortable),5)
colnames(mutsig_svsig_plot) <- c("Mutsig","cor","SVsig","p.value")

mutsig_svsig_plot[mutsig_svsig_plot$p.value > 0.05,"p.value"]  <- NA
mutsig_svsig_plot[is.na(mutsig_svsig_plot$p.value) == T,"cor"] <- NA
ggplot(data=mutsig_svsig_plot,aes(x=SVsig,y=Mutsig))+
  geom_count(aes(size = cor,color = p.value))+
  scale_colour_gradient(low = "#3182bd", high = "white")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey"))
