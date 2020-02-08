root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

output_directory_tsv <- file.path(root_dir,"scratch","sv-shatterseek","tsv")
if (!dir.exists(output_directory_tsv)) {
  dir.create(output_directory_tsv, recursive = TRUE)
}
output_directory <- file.path(root_dir,"scratch","sv-shatterseek")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}





## ===================== Load  Independent Specimen List =====================
independent_specimen_list <-  read.table(file.path(root_dir, "data","independent-specimens.wgs.primary-plus.tsv"), header = TRUE, sep = "\t")
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)
sample_names <- as.data.frame(bioid)
# ## ===================== Load Informations =====================
# info  <- read.delim(file.path(root_dir,"data","pbta-histologies.tsv"))
# info_histology <- t(info[info$Kids_First_Biospecimen_ID %in% bioid,c("Kids_First_Biospecimen_ID","short_histology","age_at_diagnosis_days")])
# bioid_medulloblastoma <- as.factor(info_histology["Kids_First_Biospecimen_ID",info_histology["short_histology",] == "Medulloblastoma"])
# 
# ## ===================== sample_names <- bioid_medulloblastoma =====================
# sample_names <- as.data.frame(bioid_medulloblastoma)


## ===================== Files and Path =====================
shatterseek_file <- read.csv(file.path(root_dir,"scratch","sv-shatterseek","PBTA_chromothripsis_newlink_region.csv"),header = TRUE)
vcf_file_path <- file.path(root_dir,"scratch","sv-vcf")

## ===================== ssmatrix save no-tra column and tra column =====================
ssmatrix_notra=matrix(0,nrow=as.numeric(nrow(sample_names)),ncol=105) ## build a ssmatrix to save no-tra column
colnames(ssmatrix_notra) <- c(paste0('svcat', 1:105,sep="")) ##column name
rownames(ssmatrix_notra) <- sample_names[,1] ## row name is sample name
ssmatrix_tra=matrix(0,nrow=as.numeric(nrow(sample_names)),ncol=7) ## build a ssmatrix to save no-tra column
colnames(ssmatrix_tra) <- c(paste0('svcat', 106:112,sep="")) ##column name
rownames(ssmatrix_tra) <- sample_names[,1] ## row name is sample name

## ===================== Remove Chrss SV and Generate Signature Input File =====================
for (i in 1:nrow(sample_names)) {
  sample <-as.character(sample_names[i,1])
  sample_vcf  <-  read.delim(file.path(vcf_file_path,paste0(sample,".tsv")))
  shatterseek_sample <-  shatterseek_file[shatterseek_file$sample == as.character(sample),]
  # out is a vcf file without chrss sv
  out <-  data.frame()
  for (j in 1:nrow(sample_vcf)) {
    chrom1 <- as.character(sample_vcf$chrom1[j])
    chrom2 <- as.character(sample_vcf$chrom2[j])
    pos1 <-  as.character(sample_vcf$pos1[j])
    pos2 <-  as.character(sample_vcf$pos2[j])
    if ((chrom1 %in% shatterseek_sample[, 1]) == TRUE &
        (chrom2 %in% shatterseek_sample[, 1]) == TRUE) {
      if (pos1 >= shatterseek_sample[shatterseek_sample$chr == chrom1, 2] &
          pos1 <= shatterseek_sample[shatterseek_sample$chr == chrom1, 3] &
          pos2 >= shatterseek_sample[shatterseek_sample$chr == chrom2, 2] &
          pos2 <= shatterseek_sample[shatterseek_sample$chr == chrom2, 3]) {
      } else {
        out <- rbind(out,sample_vcf[j,])
      }
    } else {
      out <- rbind(out,sample_vcf[j,])
    }
  }
  
  write.table(out,file.path(output_directory_tsv,paste0(sample,"_removechrss.tsv")),sep = "\t",quote = FALSE,row.names = TRUE,col.names = NA)
  
  vcf <- out
  vcf$homolen[is.na(vcf$homolen)]  <- 0
  vcf$inserlen[is.na(vcf$inserlen)]  <- 0
  ## ================== NO TRA MATRIX ========================
  vcf_notra <- vcf[vcf$SVtype !=  "TRA",]
  i2 = sample
  vcf_notra$zoomsize <- abs(vcf_notra$pos1 - vcf_notra$pos2)
  if (nrow(vcf_notra) != 0) {
    for (j in 1:nrow(vcf_notra)) {
      type <- as.character(vcf_notra$SVtype[j])
      if (type == "DEL") {
        svn = 1
      } else if (type == "DUP") {
        svn = 2
      } else if (type == "h2hINV") {
        svn = 3
      } else if (type == "t2tINV") {
        svn = 3
      }
      
      size  <- 0
      zoomsize <- vcf_notra$zoomsize[j]
      if ((zoomsize > 1) & (zoomsize <= 1000)) {
        size = 1
      } else if ((zoomsize > 1000) & (zoomsize <= 10000)) {
        size = 2
      } else if ((zoomsize > 10000) & (zoomsize <= 100000)) {
        size = 3
      } else if ((zoomsize > 100000) & (zoomsize <= 1000000)) {
        size = 4
      } else if (zoomsize > 1000000) {
        size = 5
      }
      
      homseqn <- 0
      homseq <- vcf_notra$homolen[j]
      if (homseq == 0 | homseq == 1) {
        homseqn = 1
      } else if ((homseq > 1) & (homseq <= 3)) {
        homseqn = 2
      } else if ((homseq > 3) & (homseq <= 5)) {
        homseqn = 3
      } else if (homseq > 5) {
        homseqn = 4
      }
      insl = vcf_notra$inserlen[j]
      if ((insl > 0 & insl <= 5)) {
        homseqn = 5
      } else if ((insl > 5 & insl <= 10)) {
        homseqn = 6
      } else if ((insl > 10)) {
        homseqn = 7
      }
      
      maxcell <- (svn - 1) * 35 + (size - 1) * 7 + homseqn
      ssmatrix_notra[as.character(i2), maxcell] <- as.numeric(ssmatrix_notra[as.character(i2), maxcell]) + 1
      
    }
  }
  
  
  ## ================== TRA MATRIX ========================
  vcf_tra = vcf[vcf$SVtype ==  "TRA",]
  if (nrow(vcf_tra) != 0) {
    for (k in 1:nrow(vcf_tra)) {
      homseq = vcf_tra$homolen[k]
      homseqn = 0
      if (homseq == 0 | homseq == 1) {
        homseqn = 1
      } else if ((homseq > 1) & (homseq <= 3)) {
        homseqn = 2
      } else if ((homseq > 3) & (homseq <= 5)) {
        homseqn = 3
      } else if (homseq > 5) {
        homseqn = 4
      }
      insl = vcf_tra$inserlen[k]
      if ((insl > 0 & insl <= 5)) {
        homseqn = 5
      } else if ((insl > 5 & insl <= 10)) {
        homseqn = 6
      } else if ((insl > 10)) {
        homseqn = 7
      }
      maxcell = homseqn
      ssmatrix_tra[as.character(i2), maxcell] = as.numeric(ssmatrix_tra[as.character(i2), maxcell]) + 1
    }
  }
}

ssmatrix=merge(ssmatrix_notra,ssmatrix_tra,by.x=("row.names"),by.y=("row.names"))
colnames(ssmatrix)[1]  <- "SampleName"

outputname <-  file.path(output_directory,"PBTA_SV_112matrix_withoutChrss.csv")
write.csv(ssmatrix,outputname,row.names = FALSE)
