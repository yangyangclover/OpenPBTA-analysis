
library(stringr) 
library(rprojroot)
library(readr)
library(dplyr)
library(ggplot2)

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

## =====================  Create A Subdirectory to Hold All The Output Files ===================== 
output_directory <- file.path(root_dir, "scratch","sv-vcf")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data","independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load SV File =====================
# read sv file
# filter independent-specimens
# filter "PASS" samples
sv_analysis <- read_tsv(file.path(root_dir, "data","pbta-sv-manta.tsv.gz")) %>%
  filter((Kids.First.Biospecimen.ID.Tumor %in% bioid) & (FILTER == "PASS"))

# creat small_del
keeps <- c("SV.chrom","SV.start","SV.end","SV.length","SV.type","INFO","FORMAT","Normal","Tumor","GCcontent_left","GCcontent_right","Repeats_coord_left" ,"Repeats_type_left","Repeats_coord_right","Repeats_type_right","Kids.First.Biospecimen.ID.Tumor")
small_del <- sv_analysis[,keeps]
small_del <- small_del[small_del$SV.type == "DEL",]
# add homo number
small_del$Homo <- stringr::str_match(string = small_del$INFO, pattern = "HOMLEN=[:digit:]*")
small_del$Homo <- as.numeric(gsub("HOMLEN=","",small_del$Homo))
# get support reads number
small_del_PR <-small_del[small_del$FORMAT  == "PR",]
small_del_PRSR<-small_del[small_del$FORMAT  == "PR:SR",]
small_del_PR$SUPREAD <- gsub(",","",str_match(pattern = ",[:digit:]*",small_del_PR$Tumor))
small_del_PRSR$PR <- unlist(strsplit(as.character(small_del_PRSR$Tumor),":"))[seq(1,2*nrow(small_del_PRSR),2)]
small_del_PRSR$SR <- unlist(strsplit(as.character(small_del_PRSR$Tumor),":"))[seq(2,2*nrow(small_del_PRSR),2)]
small_del_PRSR$SUPREAD <- as.numeric(gsub(",","",str_match(pattern = ",[:digit:]*",small_del_PRSR$PR))) + 
  as.numeric(gsub(",","",str_match(pattern = ",[:digit:]*",small_del_PRSR$SR)))
keeps <- colnames(small_del_PR)
small_del <- rbind(small_del_PR[,keeps],small_del_PRSR[,keeps])
# is.na >>>> 0
small_del[is.na(small_del$Homo),"Homo"] <-0

#plot
small_del$SUPREAD <- as.numeric(small_del$SUPREAD)
ggplot(data = small_del,aes(x=SV.length,y=Homo)) + 
  geom_point(aes(colour=SUPREAD),alpha=0.5,size=1.5) +
  scale_colour_gradient2(low = "gray99", high = "gray0",mid="gray80", midpoint = 10,limits=c(1,50))+
  scale_x_log10()+
  scale_y_sqrt()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 




# others
a_large=as.data.frame(small_del[small_del$SV.length > 1000 & small_del$Homo <=5,])
a_large$DEL_lENGTH_HOMO <- "length>1000_and_homo<=5"
a_large_homolarge=as.data.frame(small_del[small_del$SV.length > 1000 & small_del$Homo >5,])
a_large_homolarge$DEL_lENGTH_HOMO <- "length>1000_and_homo>5"
a_small=as.data.frame(small_del[small_del$SV.length <= 1000 & small_del$Homo <=5,])
a_small$DEL_lENGTH_HOMO <- "length<=1000_and_homo<=5"
a_small_homolarge <- as.data.frame(small_del[small_del$SV.length <= 1000 & small_del$Homo >5,])
a_small_homolarge$DEL_lENGTH_HOMO <-  "length<=1000_and_homo>5"



a <- rbind(a_large,a_large_homolarge,a_small,a_small_homolarge)


ggplot(data = a,aes(x=DEL_lENGTH_HOMO,y=SUPREAD)) + 
  geom_violin(aes(fill=factor(DEL_lENGTH_HOMO)))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 
ggplot(data = a,aes(x=DEL_lENGTH_HOMO,y=SUPREAD)) + 
  geom_boxplot(aes(fill=factor(DEL_lENGTH_HOMO)))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 


b=data.frame(names(table(a_small_homolarge$Kids.First.Biospecimen.ID.Tumor)),
  table(a_small_homolarge$Kids.First.Biospecimen.ID.Tumor))
patientgroup1 <- b[as.numeric(b$Freq)<=2,]
patientgroup1 <- factor(patientgroup1$Var1)
patientgroup2 <- b[as.numeric(b$Freq)<=4 & as.numeric(b$Freq) >2,]
patientgroup2 <- factor(patientgroup2$Var1)
patientgroup3 <- b[as.numeric(b$Freq)<=9 & as.numeric(b$Freq) >4,]
patientgroup3 <- factor(patientgroup3$Var1)
patientgroup4 <- b[as.numeric(b$Freq)<=30 & as.numeric(b$Freq) >9,]
patientgroup4 <- factor(patientgroup4$Var1)
patientgroup5 <- b[as.numeric(b$Freq)<=50 & as.numeric(b$Freq) >30,]
patientgroup5 <- factor(patientgroup5$Var1)
patientgroup6 <- b[as.numeric(b$Freq)>50,]
patientgroup6 <- factor(patientgroup6$Var1)

for (i in 1:nrow(a)) {
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup1){
    a[i,"patientgroup"] <- "patientgroup1"
  }
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup2){
    a[i,"patientgroup"] <- "patientgroup2"
  }
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup3){
    a[i,"patientgroup"] <- "patientgroup3"
  }
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup4){
    a[i,"patientgroup"] <- "patientgroup4"
  }
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup5){
    a[i,"patientgroup"] <- "patientgroup5"
  }
  if  (a[i,"Kids.First.Biospecimen.ID.Tumor"] %in% patientgroup6){
    a[i,"patientgroup"] <- "patientgroup6"
  }
}
a[is.na(a$patientgroup),"patientgroup"]<-"patientgroup1"
c <- a[a$SV.length<=1000 & a$Homo >5,]
ggplot(data = c,aes(x=patientgroup,y=SUPREAD)) + 
  geom_violin(aes(fill=factor(patientgroup)))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 
ggplot(data = c,aes(x=patientgroup,y=SUPREAD)) + 
  geom_boxplot(aes(fill=factor(patientgroup)))+
  scale_y_log10()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank())

ggplot(data = a,aes(x=SV.length,y=Homo)) + 
  geom_point(aes(colour=patientgroup),alpha=0.5,size=1.5)+
  scale_x_log10()+
  scale_y_sqrt()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 


## add disease
info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
info <- info[info$Kids_First_Biospecimen_ID %in% bioid,c("Kids_First_Biospecimen_ID","short_histology")]
a <- merge(a,info,by.x="Kids.First.Biospecimen.ID.Tumor",by.y = "Kids_First_Biospecimen_ID",all.x = T)
a$short_histology_simple <- a$short_histology
a[as.character(a[,"short_histology_simple"]) %in% names(table(a[,"short_histology"])[order(table(a[,"short_histology"]),decreasing=T)[1:10]]) ==  FALSE,"short_histology_simple"] <- "Other"

ggplot(data = a,aes(x=SV.length,y=Homo)) + 
  geom_point(aes(colour=short_histology_simple),alpha=0.5,size=1.5)+
  scale_x_log10()+
  scale_y_sqrt()+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank()) 
