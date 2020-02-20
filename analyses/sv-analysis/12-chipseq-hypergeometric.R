
setwd("~/GitHub/OpenPBTA-analysis")


library(readr)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38.masked)


# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))


genome <- BSgenome.Hsapiens.UCSC.hg38.masked

N <- maskedwidth(masks(gaps(genome$chr1))) + maskedwidth(masks(gaps(genome$chr2))) + maskedwidth(masks(gaps(genome$chr3))) +
  maskedwidth(masks(gaps(genome$chr4))) + maskedwidth(masks(gaps(genome$chr5))) + maskedwidth(masks(gaps(genome$chr6))) +
  maskedwidth(masks(gaps(genome$chr7))) + maskedwidth(masks(gaps(genome$chr8))) + maskedwidth(masks(gaps(genome$chr9))) +
  maskedwidth(masks(gaps(genome$chr10))) + maskedwidth(masks(gaps(genome$chr11))) + maskedwidth(masks(gaps(genome$chr12))) +
  maskedwidth(masks(gaps(genome$chr13))) + maskedwidth(masks(gaps(genome$chr14))) + maskedwidth(masks(gaps(genome$chr15))) +
  maskedwidth(masks(gaps(genome$chr16))) + maskedwidth(masks(gaps(genome$chr17))) + maskedwidth(masks(gaps(genome$chr18))) +
  maskedwidth(masks(gaps(genome$chr19))) + maskedwidth(masks(gaps(genome$chr20))) + maskedwidth(masks(gaps(genome$chr21))) +
  maskedwidth(masks(gaps(genome$chr22))) + maskedwidth(masks(gaps(genome$chrX))) + maskedwidth(masks(gaps(genome$chrY)))
N <- 1040147420 + 728870301 + 540864157 + 367200767 + 260556468

allfiles <- list.files(file.path(root_dir,"scratch","chipseq"),pattern = "*.gz")
allsv <- list.files(file.path(root_dir,"scratch","sv-largeall"),pattern = "*.tsv")

# function to see if breakpoint is markeable
maskable <- function(breakpoints) {
  chromnumber <- breakpoints[1]
  if(breakpoints[1] == "X") {
    chromnumber <- 23
  }
  if(breakpoints[1] == "Y") {
    chromnumber <- 24
  }
  chromnumber <- as.numeric(chromnumber)
  posnumber <- as.numeric(gsub(" ","",breakpoints[2]))
  out <- as.character(unmasked(genome[[chromnumber]])[posnumber])
  return(out)
}

# function to see if breakpoint is in chipseq
chipseq <- function(breakpoints){
  if(sum((as.character(breakpoints[1]) == bedfile[1]) & (as.numeric(breakpoints[2]) >= bedfile[2]) & (as.numeric(breakpoints[2]) <= bedfile[3]))>0){
    out <- "chipseq"
  } else {
    out <- "outchipseq"
  }
  return(out)
}


resultout <- data.frame(matrix(NA,ncol = 7,nrow = 0))
colnames(resultout) <- c("bedfile",
                         "samplesv",
                         "N",
                         "M",
                         "n",
                         "k",
                         "phyper")
count1<-0
for (filei in allfiles[5]) {
  count1 <- count1+1
  print(paste0("bedfile_number:",count1,"/386"))
  
  filepath <- file.path(root_dir,"scratch", "chipseq", filei)
  bedfile <- read_tsv(filepath, col_names=F)
  bedfile$X1 <- gsub("chr","",bedfile$X1)
  bedfile <- bedfile[bedfile$X1 %in% c(seq(1,22,1),"X","Y"),]
  
  # calculate M
  M <- 0
  for (i in 1:nrow(bedfile)) {
    chr <- paste0("chr",bedfile[i,1])
    start <- as.numeric(bedfile[i,2])
    end <- as.numeric(bedfile[i,3])
    length <- end - start + 1
    M <- M + length - str_count(unmasked(genome[[chr]])[start:end], c("N"))
  }
  
  # calculate n and  k
  count2<-0
  for (svi in allsv) {
    count2 <- count2 + 1
    print(paste0("largesvfile_number:",count2,"/797"))
    print(paste0("process:",count2+767*(count1-1),"/307642=",(count2+767*(count1-1))/307642 ))
    samplesv <- read.delim(file.path(root_dir,"scratch","sv-largeall",svi))
    samplesv <- samplesv[samplesv$largeall == "largeall",]
    if (nrow(samplesv) > 0){
      breakpoints1 <- data.frame(samplesv$chrom1,samplesv$pos1)
      breakpoints2 <- data.frame(samplesv$chrom2,samplesv$pos2)
      breakpoints <- rbind(breakpoints1, setNames(breakpoints2, names(breakpoints1)))
      breakpoints <- breakpoints[breakpoints$samplesv.chrom1 != "M",]
      breakpoints$maskable <- apply(breakpoints,1,maskable)
      breakpoints$chipseq <- apply(breakpoints,1,chipseq)
      
      n <- nrow(breakpoints[breakpoints$maskable %in% c("A","T","C","G"),])
      k <- nrow(breakpoints[breakpoints$chipseq == "chipseq",])
      
      # phyper
      phyper_value <- 1-phyper(k-1,M,N-M,n)
      
      # write phyper_value <  0.05
      if (phyper_value >0) {
        resultout <- rbind(resultout,data.frame(
          "bedfile" = filei,
          "samplesv" = svi,
          "N"=N,
          "M"=M,
          "n"=n,
          "k"=k,
          "phyper"=phyper_value
        ))
      }
    }
  }
}
write.csv(resultout,file.path(root_dir,"scratch","chipseq","phyper_chipseq_sample_result(SK-N-SH_H3K27ac).csv"),quote = F)
