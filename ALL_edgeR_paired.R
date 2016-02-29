rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_ALL/HtSeq_reverse/")

library("edgeR")
library("matrixStats")

###################
# read HtSeq data #
###################

# samples, tech reps, unique counts
all.files <- c(
  "ALL_SORT_B1T1_HtseqCounts.txt",
  "ALL_SORT_B1T2_HtseqCounts.txt",
  "ALL_SORT_B2T1_HtseqCounts.txt",
  "ALL_SORT_B2T2_HtseqCounts.txt",
  "ALL_SORT_B3T1_HtseqCounts.txt",
  "ALL_SORT_B3T2_HtseqCounts.txt",
  "ALL_SORT_B4T1_HtseqCounts.txt",
  "ALL_SORT_B4T2_HtseqCounts.txt",
  "NO_SORT_B1T1_HtseqCounts.txt",
  "NO_SORT_B1T2_HtseqCounts.txt",
  "NO_SORT_B2T1_HtseqCounts.txt",
  "NO_SORT_B2T2_HtseqCounts.txt",
  "NO_SORT_B3T1_HtseqCounts.txt",
  "NO_SORT_B3T2_HtseqCounts.txt",
  "NO_SORT_B4T1_HtseqCounts.txt",
  "NO_SORT_B4T2_HtseqCounts.txt",
  "YFP_NEG_B1T1_HtseqCounts.txt",
  "YFP_NEG_B1T2_HtseqCounts.txt",
  "YFP_NEG_B2T1_HtseqCounts.txt",
  "YFP_NEG_B2T2_HtseqCounts.txt",
  "YFP_NEG_B3T1_HtseqCounts.txt",
  "YFP_NEG_B3T2_HtseqCounts.txt",
  "YFP_NEG_B4T1_HtseqCounts.txt",
  "YFP_NEG_B4T2_HtseqCounts.txt",
  "YFP_POS_B1T1_HtseqCounts.txt",
  "YFP_POS_B1T2_HtseqCounts.txt",
  "YFP_POS_B2T1_HtseqCounts.txt",
  "YFP_POS_B2T2_HtseqCounts.txt",
  "YFP_POS_B3T1_HtseqCounts.txt",
  "YFP_POS_B3T2_HtseqCounts.txt"
)

# reads mapped to reverse strand


# read files and merge into data frame all.counts
c <- 1
for (filename in all.files) {
  if(c == 1){ # if it is the first file just read file
    all.counts <- read.table(filename,sep="\t",stringsAsFactors=F)
    names(all.counts) <- c("ID", c)
  }
  else{ # else merge the other files
    tmp <- read.table(filename,sep="\t",stringsAsFactors=F)
    names(tmp) <- c("ID",c)
    all.counts <- merge(all.counts,tmp,by="ID")
  }
  cat(filename, "read.\n")
  c <- c+1
}
# store results
all.counts.storage <- all.counts

# remove rows that do not correspond to Tair loci (HtSeq stats)
all.counts <- all.counts.storage[-c(1,2,33605,33606,33607),]

# merge tech replicates
merged.counts <- data.frame(all.counts[,1],stringsAsFactors=F)
for (i in 1:15) {
  merged.counts[,i+1] <- all.counts[,2*i]+all.counts[,1+2*i]
}
all.counts <- merged.counts

# set variable names
names(all.counts) <- c("GeneID", "ALL_SORT_B1T1", "ALL_SORT_B2T1", "ALL_SORT_B3T1", "ALL_SORT_B4T1", "NO_SORT_B1T1", "NO_SORT_B2T1", "NO_SORT_B3T1", "NO_SORT_B4T1", "YFP_NEG_B1T1", "YFP_NEG_B2T1", "YFP_NEG_B3T1", "YFP_NEG_B4T1", "YFP_POS_B1T1", "YFP_POS_B2T1", "YFP_POS_B3T1")
str(all.counts)

# use gene length (likely derived from major transcript)
len <- read.table("tair10.whole.genelength.txt",sep="\t")
names(len) <- c("GeneID","length")
all.counts <- merge(all.counts,len,by=c(1))
# re-order variables 
all.counts <- all.counts[,c(1,17,2:16)]

# generate count table YFP_NEG, YFP_pos, drop fourth biological replicate
cm <- as.matrix(all.counts[,-c(1:10,14)])

# If you want to do paired analysis with ALL/NO_sort use this insted
# generate reduced count table:  ALL_SORT vs NO_SORT
cm <- as.matrix(all.counts[,-c(1,2,11:17)])

# set row names
rownames(cm) <- all.counts[,1]

# build DGEList
y <- DGEList(counts = cm)
# filtering low expressed genes
min.cpm <- 2
n.min.cpm <- 3
keep <- rowSums( cpm(y)>min.cpm ) >= n.min.cpm
table(keep)
y <- y[keep,]
dim(y)

# normalization
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y$samples
plotMDS(y)

# build design matrix yfp_neg vs yfp_pos
run <- factor(c(1,2,3,1,2,3)) 
tissue <- factor(c("YFP_NEG", "YFP_NEG", "YFP_NEG", "YFP_POS", "YFP_POS", "YFP_POS"))
sample <- c("YFP_NEG_B1T1", "YFP_NEG_B2T1", "YFP_NEG_B3T1", "YFP_POS_B1T1", "YFP_POS_B2T1", "YFP_POS_B3T1")
data.frame(sample, run, tissue)

# build design matrix All_Sort vs NO_sort
run <- factor(c(1,2,3,4,1,2,3,4)) 
tissue <- factor( c("ALL_SORT", "ALL_SORT", "ALL_SORT", "ALL_SORT", "NO_SORT", "NO_SORT", "NO_SORT", "NO_SORT"))  
sample <- c("ALL_SORT_B1T1", "ALL_SORT_B2T1", "ALL_SORT_B3T1", "ALL_SORT_B4T1", "NO_SORT_B1T1", "NO_SORT_B2T1", "NO_SORT_B3T1", "NO_SORT_B4T1")  
data.frame(sample, run, tissue)
 

design <- model.matrix(~run+tissue)
rownames(design) <- colnames(y)
design

# estimate dispersion
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
etable <- topTags(lrt, n=nrow(lrt$table), adjust.method="BH")
#etable <- etable$table[etable$table$FDR<fdr.t,]
etable <- etable$table[etable$table$FDR<=1,]
write.table( etable[,], file=paste("YFP_NEG_POS_paired", min.cpm,n.min.cpm,sep="."), row.names=TRUE)


