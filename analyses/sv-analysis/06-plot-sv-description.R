library(reshape2)
library(ggplot2)
library(ggthemes)



# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

## =====================  Create A Subdirectory to Hold All The Output Files ===================== 
output_directory <- file.path(root_dir, "scratch","sv-description","plots")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data/independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load tsv file and generate sv_matrix ====================
tsv_file_path <- file.path(root_dir,"scratch","sv-vcf")
# tsv_file_path <- file.path(root_dir,"scratch","sv-shatterseek","tsv")

sv_df <-  data.frame()
for (j in bioid) {
  sample <- j
  if (file.exists(file.path(tsv_file_path,paste0(sample,".tsv")))){
  # if (file.exists(file.path(tsv_file_path,paste0(sample,"_removechrss.tsv")))){
    sample_tsv  <-  read.delim(file.path(tsv_file_path,paste0(sample,".tsv")))
    # sample_tsv  <-  read.delim(file.path(tsv_file_path,paste0(sample,"_removechrss.tsv")))
    sv_df_add <- data.frame(sample = j,
                            DEL = nrow(sample_tsv[sample_tsv$SVtype == "DEL",]),
                            DUP = nrow(sample_tsv[sample_tsv$SVtype == "DUP",]),
                            INV = nrow(sample_tsv[sample_tsv$SVtype == "t2tINV",]) + nrow(sample_tsv[sample_tsv$SVtype == "h2hINV",]),
                            TRA = nrow(sample_tsv[sample_tsv$SVtype == "TRA",]))
    sv_df <- rbind(sv_df,sv_df_add)
  }
} 
sv_df$sum <- rowSums(sv_df[,c("DEL","DUP","INV","TRA")])

## ===================== Plot sv number =====================
order <- order(sv_df$sum,decreasing = T)
sv_df <- sv_df[order,]
sv_df_plot <-  melt(sv_df,id =  c("sample"), measure.vars = c("DEL","DUP","INV","TRA"),value.name = "Number",variable.name = "SVtype")
sv_df_plot$sample <- factor(sv_df_plot$sample,levels = unique(sv_df$sample))
title <- "SV Number_PBTA"
ggplot(sv_df_plot,
       aes(sample,Number,fill = SVtype)) + 
  geom_bar(stat = "identity",position="stack",width = 1) + 
  ggtitle(title) +
  labs(y = "Number",  x = "") +
  scale_fill_wsj("colors6", "")  +
  theme(axis.text.x = element_blank(),
        panel.grid.major= element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(title = NULL))
ggsave(file.path(output_directory,paste0(title,".pdf")),plot = last_plot(), device = NULL)

## ===================== Plot sv percent =====================
order <- order(sv_df$TRA/sv_df$sum,sv_df$DEL/sv_df$sum,sv_df$DUP/sv_df$sum,sv_df$INV/sv_df$sum,decreasing=T)
sv_df <- sv_df[order,]
sv_df_plot <-  melt(sv_df,id =  c("sample"), measure.vars = c("DEL","DUP","INV","TRA"),value.name = "Number",variable.name = "SVtype")
sv_df_plot$sample <- factor(sv_df_plot$sample,levels = unique(sv_df$sample))
title <- "SV Percent_PBTA"
ggplot(sv_df_plot,
       aes(sample,Number,fill = SVtype)) + 
  geom_bar(stat = "identity",position="fill",width = 1) + 
  ggtitle(title) +
  labs(y = "Percent",  x = "") +
  scale_fill_wsj("colors6", "")  +
  theme(axis.text.x = element_blank(),
        panel.grid.major= element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(title = NULL))
ggsave(file.path(output_directory,paste0(title,".pdf")),plot = last_plot(), device = NULL)

## ===================== Plot sv-info =====================
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
sv_info <- merge(sv_df,info[,c("Kids_First_Biospecimen_ID","short_histology")],by.x="sample",by.y = "Kids_First_Biospecimen_ID")
order <- order(sv_info$sum,decreasing = T)
sv_info <- sv_info[order,]
sv_info$sample <- factor(sv_info$sample,levels = sv_info$sample)

ggplot(sv_info,aes(x=sample,y=sum,color=short_histology))+
  geom_point(size=1)+
  facet_grid(. ~ short_histology,scales="free_x")+
  scale_y_log10()+
  scale_fill_manual("legend_title") +
  labs(col="Tumor Type")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=10), # cancel x text
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=10),
        axis.ticks.length=unit(0,"cm"), # cancel ticks
        panel.background = element_rect(fill = "gray94",colour = "gray94",size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray94"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray94"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=9),
        legend.position="bottom",
        axis.line=element_line(colour="black"), # cancel legend, change x & y axis lines to black
        axis.line.x=element_line(size=2),
        axis.line.y=element_line(size=1),
        line=element_line(size=1.5)) #line size of x is 2, line size of y is 1, scale line size is 1.5
ggsave(file.path(output_directory,"sv_number_across_tumors.pdf"),dpi = 2000,width = 35, height = 20, units = "cm")
