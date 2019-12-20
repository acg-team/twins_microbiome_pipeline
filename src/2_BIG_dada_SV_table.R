# @Alex
############# Step 2 #################
# Create a joint SV table for all samples together (analogous to OTU table)
# 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# https://github.com/benjjneb/dada2/tree/master/R


##### init: load packages and set path  
getwd()
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")
packageVersion("dada2")

### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file)) 




###### 1: Files preparation / Quality check ########################

# get all sample's file names to be processed
# either by scanning the data folder
fns <- sort(list.files(data_path, full.names = TRUE))
print(paste("Samples (reverse and forward) left BEFORE filtering::", length(fns)))

# keep only those families who have 4 samples 
# NOTE!!!:: REMOVE it if you need full file list!!! only 700 samples are left now!
samples_to_be_kept <- df.metadata.4timepoints$file
fns <- fns[grep(paste(samples_to_be_kept, collapse="|"), fns)]
print(paste("Samples (reverse and forward) left after filtering::", length(fns)))

# split all filenames into Forward and Reverse reads files
fnFs <- fns[grepl("_1.fastq.gz", fns)]
fnRs <- fns[grepl("_2.fastq.gz", fns)]

# retrieve the sample names from file names and keep them in a list ("ERR1382288" etc)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

# check if there are both forward and reverse files
for (idx in seq_along(fnFs)){
  if(sample.names[idx]!=sample.namesR[idx]) 
    {
    print(idx)
    print(fnFs[[idx]])
    print(fnRs[[idx]])
    stop("Files: Forward and reverse files do not match!")
  }
}


### QUALITY: ##############
# check FastQC plots for quality 
# some of them are in bad quality #3
ii <- seq(from=1,to=5,by=1)  #length(fnFs)
for(i in ii) {
 print(dada2::plotQualityProfile(fnFs[i]) + ggtitle(paste("Fwd:", sample.names[i])))
 print(dada2::plotQualityProfile(fnRs[i]) + ggtitle(paste("Rev:", sample.names[i])))
}


### FILTERING  ##########################
# trim and put into filtered folder
# https://github.com/benjjneb/dada2/tree/master/R
QUALITY_THRESHOLD <- 20  # Phred

filtFs <- file.path(filt_path, basename(fnFs))  # names for filtered forwards reads
filtRs <- file.path(filt_path, basename(fnRs))  # names for filtered reverse reads

names(filtFs) <- sample.names
names(filtRs) <- sample.names

if(length(fnFs) != length(fnRs)) stop("BEFORE: Forward and reverse files do not match!")

tic()
# quality filtering and trimming
# TODO: need an expert advice! Or experiment with eliminating bad quality reads !
# need to deside on: 
# - trimLeft/truncLen; 
# - maxEE - maximum number of “expected errors” allowed in a read
# - truncQ;
# https://academic.oup.com/bioinformatics/article/31/21/3476/194979
for(i in seq_along(fnFs)) {
  print(i)
  print(fnFs[[i]])
  print(fnRs[[i]])
  dada2::filterAndTrim( fwd=fnFs[[i]], filt=filtFs[[i]],
                        rev=fnRs[[i]], filt.rev=filtRs[[i]],
                        #trimLeft=c(3,3), truncLen=c(247,235), # warning: No reads passed the filter ?!
                        maxEE=c(2,5), maxN=0, truncQ=QUALITY_THRESHOLD,  rm.phix=TRUE,
                        compress=TRUE, verbose=TRUE, multithread=TRUE
  )
}
toc()  # 5839sec = 1.5 hours


# check quality afterward - if nessesary
ii <- seq(from=1,to=5,by=1)  #length(fnFs)
for(i in ii) { 
  print(plotQualityProfile(fnFs[i]) + ggtitle(paste("Fwd:", sample.names[i]))) 
  print(plotQualityProfile(filtFs[i]) + ggtitle(paste("Fwd_After:", sample.names[i]))) 
  
  print(plotQualityProfile(fnRs[i]) + ggtitle(paste("Rev:", sample.names[i]))) 
  print(plotQualityProfile(filtRs[i]) + ggtitle(paste("Rev_After:", sample.names[i]))) 
}


########################################################################
###### BIG DATA sequential workflow in linear time #####################
if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match!")


### LEARN ERROR RATES ###############
# https://github.com/benjjneb/dada2/issues/155
# for large data sets error rates should be estimated on subset of data - change to 40, here it is ok
tic()
errF <- learnErrors(filtFs, nreads=2e6, multithread = TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nreads=2e6, multithread = TRUE, randomize=TRUE)
print("Error rate calculation time:")
toc() # 9806.114 sec / 2.7h = now 2509sec

## plot error rates for control
plotErrors(errF)
plotErrors(errR)
save(errF, errR, file=file.path(files_intermediate, dada.err.file)) 


### Big Data dada2 pipeline
# do it in a loop because we have a huge dataset (BigData workflow)
mergers <- list()
counter <- 0

tic()
for (sam in sample.names) {
  tic()
  ### DEREPLICATION 
  derepF <- derepFastq(filtFs[[sam]])
  derepR <- derepFastq(filtRs[[sam]])
  
  ### SAMPLE INFERENCE 
  dadaF <- dada(derepF, err=errF, pool=TRUE, multithread = TRUE)
  dadaR <- dada(derepR, err=errR, pool=TRUE, multithread = TRUE)
  
  merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  mergers[[sam]] <- merger
  toc()
  
  counter <- counter+1
  cat(counter, "...", sam, " Done.\n")
  cat("---------- \n")
  save(mergers, file=file.path(files_intermediate, mergers.file)) 
}
print("Total time of sample inference:")
toc() # 8669sec



### SEQTAB: constructs a sequence table (analogous to an OTU table) from the list of samples.
tic()
seqtab.all <- dada2::makeSequenceTable(mergers, orderBy = "abundance")
print("Total time of makeSequenceTable:")
toc() # 5sec

### remove chimeras ----
seqtab <- dada2::removeBimeraDenovo(seqtab.all, verbose = TRUE)
save(seqtab, file=file.path(files_intermediate, seqtab.file)) 


### Extract sample names and save them separatelly (for futher Python data analysis)
seqtab.samples.names = rownames(seqtab)
save(seqtab.samples.names, file=file.path(files_intermediate, seqtab.snames.file)) 


### TODO:
# sequence variants have very long names, so substitute them with short ones (Seq1, Seq2 ...) and create mathing table
