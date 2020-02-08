library(ComplexHeatmap)
library(ggplot2)
library(ggthemes)


## ===================== Root Directory =====================
# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


## =====================  Create A Subdirectory to Hold All The Output Files ===================== 
# save output files in "scratch/sv-shatterseek" and "analyses/sv-analysis/tables/sv-shatterseek"
output_directory <- file.path(root_dir, "scratch","sv-shatterseek")


## ===================== Load Shatterseek_link file ===================== 
list <- read.csv(file.path(root_dir,"scratch","sv-shatterseek","PBTA_chromothripsis_newlink_region.csv"))
chrsssampleid <- as.data.frame(unique(list[,"sample"]))
colnames(chrsssampleid) <- "sample"
chrsssampleid$chrss <- 1

## ===================== Load  Independent specimen list =====================
independent_specimen_list <- read.table(file.path(root_dir,"data","independent-specimens.wgs.primary-plus.tsv"),header = TRUE,sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load Info ========================
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))

## ===================== Load  SV ============================
tsv_file_path <- file.path(root_dir,"scratch","sv-vcf")
sv_df <-  data.frame()
for (j in bioid) {
  sample <- j
  if (file.exists(file.path(tsv_file_path,paste0(sample,".tsv")))){
    sample_tsv  <-  read.delim(file.path(tsv_file_path,paste0(sample,".tsv")))
    sv_df_add <- data.frame(sample = j,
                            DEL = nrow(sample_tsv[sample_tsv$SVtype == "DEL",]),
                            DUP = nrow(sample_tsv[sample_tsv$SVtype == "DUP",]),
                            INV = nrow(sample_tsv[sample_tsv$SVtype == "t2tINV",]) + nrow(sample_tsv[sample_tsv$SVtype == "h2hINV",]),
                            TRA = nrow(sample_tsv[sample_tsv$SVtype == "TRA",]))
    sv_df <- rbind(sv_df,sv_df_add)
  }
} 
sv_df$sum <- rowSums(sv_df[,c("DEL","DUP","INV","TRA")])

## ===================== Creat sample_chrss_info df =======================
# creat df
sample_chrss_info <- data.frame(matrix(NA,nrow = length(bioid)))
colnames(sample_chrss_info) <- "ID"
sample_chrss_info$chrss <- as.factor(sample_chrss_info$chrss)

# ID is sample id
sample_chrss_info$ID <- bioid
# add chrss info
sample_chrss_info <- merge(sample_chrss_info,chrsssampleid,by.x = "ID",by.y = "sample",all.x = T)
sample_chrss_info[is.na(sample_chrss_info$chrss),"chrss"] <- "no-chromothripsis"
sample_chrss_info[sample_chrss_info$chrss == 1,"chrss"] <- "chromothripsis"
# add other info
sample_chrss_info <-  merge(sample_chrss_info,
                            info[,c("Kids_First_Biospecimen_ID","age_at_diagnosis_days","short_histology","germline_sex_estimate")],
                            by.x = "ID",by.y = "Kids_First_Biospecimen_ID",all.x = T)
sample_chrss_info$age_at_diagnosis_days <- as.numeric(as.character(sample_chrss_info$age_at_diagnosis_days))
# add sv number
sample_chrss_info <- merge(sample_chrss_info,sv_df,by.x="ID",by.y="sample",all.x = T)
sample_chrss_info[,c("DEL","DUP","INV","TRA")] <- as.matrix(sample_chrss_info[,c("DEL","DUP","INV","TRA")] )
# add sv percent
sample_chrss_info <- cbind(sample_chrss_info,t(apply(sample_chrss_info[,c("DEL","DUP","INV","TRA")], 1, function(x) x/sum(x))))
colnames(sample_chrss_info)[11:14] <- c("DEL_per","DUP_per","INV_per","TRA_per")
# add top 10 disease
sample_chrss_info$short_histology_simple <- sample_chrss_info$short_histology
sample_chrss_info[as.character(sample_chrss_info[,"short_histology_simple"]) %in% names(table(sample_chrss_info[,"short_histology"])[order(table(sample_chrss_info[,"short_histology"]),decreasing=T)[1:10]]) ==  FALSE,"short_histology_simple"] <- "Other"
# add tumor number
tumornumber <- as.data.frame(table(sample_chrss_info[,"short_histology"])[order(table(sample_chrss_info[,"short_histology"]),decreasing=T)])
colnames(tumornumber)<-c("Tumor","TumorNum")
sample_chrss_info<- merge(sample_chrss_info,tumornumber,by.x = "short_histology",by.y="Tumor",all.x = T)
# add tumor number simple
tumornumber_simple <- as.data.frame(table(sample_chrss_info[,"short_histology_simple"])[order(table(sample_chrss_info[,"short_histology_simple"]),decreasing=T)])
colnames(tumornumber_simple)<-c("Tumor","TumorNumSimple")
sample_chrss_info<- merge(sample_chrss_info,tumornumber_simple,by.x = "short_histology_simple",by.y="Tumor",all.x = T)

## =================== Plot chrss_info======================
disease_col <- c("blue","darkorange","cornsilk","brown3","skyblue","limegreen","cadetblue","gold","deeppink","yellow")
zero_row_mat <- matrix(nrow = 0, ncol = nrow(sample_chrss_info))

column_up <- HeatmapAnnotation("Chromothripsis"=sample_chrss_info[,"chrss"],
                               "Age"=sample_chrss_info[,"age_at_diagnosis_days"],
                               "Sex"=sample_chrss_info[,"germline_sex_estimate"],
                               "Disease"=sample_chrss_info[,"short_histology_simple"],
                               col = list(Chromothripsis=c("no-chromothripsis"="white","chromothripsis"="firebrick3"),
                                          Sex=c("Female"="violet","Male"="snow2"),
                                          Disease = setNames(disease_col, unique(info_histology["short_histology_simple",]))),
                               show_legend = c(F,T,T,T))

column_down <- HeatmapAnnotation(
  sv_number = anno_barplot(sample_chrss_info[,c("DEL","DUP","INV","TRA")],
                            height = unit(5, "cm"),
                            width = 1,
                            gp = gpar(col=c("firebrick3","steelblue4","gold","forestgreen")))
)



ht <- Heatmap(zero_row_mat,
              top_annotation = column_up,
              bottom_annotation = column_down,
              column_order = order(sample_chrss_info$sum,decreasing = T),
              show_heatmap_legend = T,
              show_column_names = FALSE,
              #clustering_distance_columns = "pearson" # ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
)

draw(ht,
     show_annotation_legend = TRUE,
     annotation_legend_side = "left")

## =================== Plot chrss across tumors======================
# order by tumornum
sample_chrss_info <- sample_chrss_info[order(sample_chrss_info$TumorNum,decreasing = T),]
sample_chrss_info$short_histology <- factor(sample_chrss_info$short_histology,levels = unique(sample_chrss_info$short_histology))
# creat chrss_tumor
chrss_tumor <- as.data.frame(matrix(NA,nrow = length(unique(sample_chrss_info$short_histology))*2,ncol = 3))
colnames(chrss_tumor) <- c("Tumor","chrss","Freq")
chrss_tumor$Tumor <- rep(unique(sample_chrss_info$short_histology),2)
chrss_tumor$chrss <- rep(c("no-chromothripsis","chromothripsis"),each=length(unique(sample_chrss_info$short_histology)))
chrss_tumor$Freq <- c(table(sample_chrss_info[,c("chrss","short_histology")])[2,1:length(unique(sample_chrss_info$short_histology))],table(sample_chrss_info[,c("chrss","short_histology")])[1,1:length(unique(sample_chrss_info$short_histology))])

ggplot(data = chrss_tumor,aes(x=Tumor,y=Freq,fill=chrss)) + 
  geom_bar(stat = "identity") +
   # geom_text(aes(label = Freq),size = 3, position = position_stack(vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.grid=element_blank(),panel.border=element_blank(),   #remove grid and border
        panel.background = element_blank(), #remove background
        axis.line = element_line(colour = "black")) +
  ylab("Number of Patient") +
  xlab("Tumor Types") +
  scale_fill_manual("legend", values = c("chromothripsis" = "firebrick3", "no-chromothripsis" = "dodgerblue3"))


## =================== Plot Disease Distribution =================
# order by tumornum
tumornumber <- tumornumber[order(tumornumber$TumorNum,decreasing = F),]
tumornumber$Tumor <- factor(tumornumber$Tumor,levels = unique(tumornumber$Tumor))

ggplot(tumornumber,aes(x = Tumor, y = TumorNum)) + 
  geom_bar(stat = "identity",width = 0.8,col="firebrick3",fill="firebrick3") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 15),
        panel.grid=element_blank(),panel.border=element_blank(),   #remove grid and border
        panel.background = element_blank(), #remove background
        axis.line = element_line(colour = "black")) +
  geom_text(aes(label = TumorNum), col="black",nudge_y = 8) +
  coord_flip() 

## =================== Plot sv across tumors ===================
# order by tumornum
sample_chrss_info <- sample_chrss_info[order(sample_chrss_info$TumorNum,sample_chrss_info$sum ,decreasing = T),]
sample_chrss_info$short_histology <- factor(sample_chrss_info$short_histology,levels = unique(sample_chrss_info$short_histology))
sample_chrss_info$ID <- factor(sample_chrss_info$ID,levels = unique(sample_chrss_info$ID))
# background color
backgroundcolor <- "grey98"
# plot
ggplot(sample_chrss_info,aes(x=ID,y=sum,color=short_histology))+
  geom_point(size=1.5)+
  facet_grid(. ~ short_histology,scales="free_x")+
  scale_y_log10()+
  scale_fill_manual("legend_title") +
  xlab("Tumor Type")+
  ylab("SV Nmuber")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=10), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=10),
        axis.ticks=element_blank(), # cancel ticks
        axis.line=element_blank(),
        panel.background = element_rect(fill = backgroundcolor,colour = backgroundcolor,size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = backgroundcolor), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = backgroundcolor),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        legend.position="bottom",
        line=element_line(size=1.5)) #line size of x is 2, line size of y is 1, scale line size is 1.5
