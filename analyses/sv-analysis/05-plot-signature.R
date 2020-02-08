## =====================  Load Packages =====================
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(randomcoloR)
library(ggpubr)
## =====================  Create A Subdirectory to Hold All The Output Files ===================== 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_directory <- file.path(root_dir, "scratch","sv-signature","plots","sigprofiler")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


## =====================  Import Signature Files ===================== 
# method can be "sigprofiler"/"signeR"/"NMF"
# file name can be changed
method <- "Sigprofiler"

# load mutation_sig_file
signature_path <-  file.path("D:","Yang","Downloads","SigProfiler_2_5_1_7","SigProfiler_2_5_1_7","output","PBTA_PASS_v7.112","text")
mutation_sig_file <- "res_PBTA_PASS_v7.112_signature_patterns_for_5_sigs.csv"
mutation_sig <- read.csv(file.path(signature_path,mutation_sig_file))
mutation_sig_analysis <- data.frame(SVtype = mutation_sig$Mutation.Type,
                                    SVsubtype = paste(mutation_sig$Mutation.Type,mutation_sig$Mutation.Subtype,sep = ":"))
mutation_sig_analysis <- cbind(mutation_sig_analysis,mutation_sig[,c(3:ncol(mutation_sig))])

# load sample_sig_file
sample_sig_file <- "res_PBTA_PASS_v7.112_signature_activities_for_5_sigs.csv"
sample_sig <- read.csv(file.path(signature_path,sample_sig_file))
sample_sig_analysis <- sample_sig
sample_sig_analysis$sum <- rowSums(sample_sig_analysis[,c(2:ncol(sample_sig_analysis))])

## =================== Load info ======================
# load information
information <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
information$Sample.Names <- substr(information$Kids_First_Biospecimen_ID,4,11)

## =================== Plot Sample_sig_plot ==============
# order [1 out of 2]
sample_sig_analysis <- sample_sig_analysis[order(sample_sig_analysis$sum,decreasing = TRUE),]
orderall <- unique(sample_sig_analysis$Sample.Names)
# total number of signature
signumber <-  ncol(sample_sig) - 1
# melt sample_sig_plot
sample_sig_plot <- melt(sample_sig_analysis[,colnames(sample_sig_analysis) != "sum"],value.name="number", variable.name="sig", na.rm=TRUE)
# order sample.names as orderall
sample_sig_plot$Sample.Names <- factor(sample_sig_plot$Sample.Names,  orderall)
# add GenomeSize AND MutePerMb
sample_sig_plot$GenomeSize <- 3088286401
sample_sig_plot$MutePerMb <- (sample_sig_plot$number/sample_sig_plot$GenomeSize) * (10^6)
# add some info, like short_histology
sample_sig_plot <- merge(sample_sig_plot,information[,colnames(information) %in% c("short_histology",
                                                                                   "Sample.Names")],
                         by.x = "Sample.Names",by.y = "Sample.Names")
# functions [plot_sig]
plot_sig <- function(df,label,position) {
  wide = 1
  if (nrow(df)/signumber < 100){
    wide  = 0.8
  }
  if(position == "stack") {
    ylabname = paste0("Value of Signature by ",method)}
  if (position == "fill") {
    ylabname = paste0("Percent of Signature by ",method)}
  main <- ggplot(data=df, aes(x=Sample.Names, y=number, fill=sig, width= wide)) +
    geom_bar(stat="identity",position=position,width= wide)+
    scale_fill_manual(values=c("Signature.1" = "firebrick3", 
                               "Signature.2" = "dodgerblue3",
                               "Signature.3" = "coral",
                               "Signature.4" = "darkorchid",
                               "Signature.5" = "darkgreen")) +
    ylab(ylabname) +
    xlab("Samples") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = rel(2)),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.ticks = element_blank(),
      panel.grid=element_blank(),panel.border=element_blank(),   #remove grid and border
      panel.background = element_blank(), #remove background
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = "top"
    )
  filename  <-  paste0("signatures_",position,".pdf")
  ggsave(file.path(output_directory,paste(label,filename,sep = "_")),width = 50, height = 30,units = "cm",dpi = 1500)
}
# Plot All in plot_sig_number
plot_sig(sample_sig_plot,label = "All",position = "stack")
# Plot Sub in plot_sig_number
for (i in unique(sample_sig_plot$short_histology)) {
  subdf <- sample_sig_plot[sample_sig_plot$short_histology == i,]
  plot_sig(subdf,label = i,position = "stack")
}
# order
sample_sig_analysis <- sample_sig_analysis[order(sample_sig_analysis$Signature.2/sample_sig_analysis$sum,
                                                 sample_sig_analysis$Signature.4/sample_sig_analysis$sum,
                                                 sample_sig_analysis$Signature.5/sample_sig_analysis$sum,
                                                 sample_sig_analysis$Signature.1/sample_sig_analysis$sum,
                                                 sample_sig_analysis$Signature.3/sample_sig_analysis$sum,
                                                 decreasing = TRUE),]
orderall <- unique(sample_sig_analysis$Sample.Names)
# Plot All in plot_sig_percent
plot_sig(sample_sig_plot,label = "All",position = "fill")
# Plot Sub in plot_sig_percent
for (i in unique(sample_sig_plot$short_histology)) {
  subdf <- sample_sig_plot[sample_sig_plot$short_histology == i,]
  plot_sig(subdf,label = i,position = "fill")
}


## =====================  Plot SV Mutation Signatures ===================== 
cols = brewer.pal(signumber, "Set1")
# melt mutation_sig_plot
mutation_sig_plot <- melt(mutation_sig_analysis,value.name="percent", variable.name="sig", na.rm=TRUE)
# order SVtype adn SVsubtype
mutation_sig_plot$SVtype <- factor(mutation_sig_plot$SVtype,levels = unique(mutation_sig_plot$SVtype))
mutation_sig_plot$SVsubtype <- factor(mutation_sig_plot$SVsubtype,levels = unique(mutation_sig_plot$SVsubtype))
# plot
ggplot(data = mutation_sig_plot, aes(x=SVsubtype,y=percent,fill=SVtype)) + 
  geom_col() +
  ylab("Percent")+
  scale_fill_manual(values=c("DEL" = "firebrick3", 
                             "DUP" = "dodgerblue3",
                             "INV" = "darkgreen",
                             "TRA" = "darkorchid")) + 
  facet_wrap("sig", ncol = 1, scales = "free") +
  theme(axis.text.x = element_text(angle = 90,size = 7 ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.background = element_rect(color = "grey99"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),) 
ggsave(file.path(output_directory,"mutation_sig_plot_figure.pdf"),width = 50, height = 30,units = "cm",dpi = 1500)

