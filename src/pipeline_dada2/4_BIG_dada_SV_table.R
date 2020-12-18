# @Alex
############# Step 2 #################
# Create a joint SV table for all samples together (analogous to OTU table)
# 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# https://github.com/benjjneb/dada2/tree/master/R


##### init: load packages and set path  
packageVersion("dada2")

### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file)) 


### FILTERING  ##########################
# trim and put into filtered folder
# https://github.com/benjjneb/dada2/tree/master/R
print(paste("--------   Filtering parameters : ", dada_param, " ---------" ))

# read parameters
QUALITY_THRESHOLD <- dada_param$QUALITY_THRESHOLD #18  # Phred
maxEE <- dada_param$maxEE   # c(5,5)
trimLeft <- dada_param$trimLeft
trimRight <- dada_param$trimRight
truncLen <- dada_param$truncLen

folder.suffix <- paste0(
  conf$dataset, "_", conf$pipeline, 
  "_Q", dada_param$QUALITY_THRESHOLD, 
  "_mEE", dada_param$maxEE[1], dada_param$maxEE[2], 
  "_trL", dada_param$trimLeft[1], dada_param$trimLeft[2],
  "_trR", dada_param$trimRight[1], dada_param$trimRight[2],
  "_truncLn", dada_param$truncLen[1], "_", dada_param$truncLen[2]
)

print(folder.suffix)

# form a suffics for filtered reads names
filt_path.suf <- paste0(filt_path, "/", folder.suffix)
  
filtFs <- file.path(filt_path.suf, basename(fnFs))  # names for filtered forwards reads
filtRs <- file.path(filt_path.suf, basename(fnRs))  # names for filtered reverse reads

names(filtFs) <- sample.names
names(filtRs) <- sample.names

if(length(fnFs) != length(fnRs)) stop("BEFORE: Forward and reverse files do not match!")

# quality filtering and trimming
# TODO: need an expert advice! Or experiment with eliminating bad quality reads !
# need to deside on: 
# - trimLeft/truncLen; 
# - maxEE - maximum number of “expected errors” allowed in a read
# - truncQ;
# https://academic.oup.com/bioinformatics/article/31/21/3476/194979
# https://github.com/benjjneb/dada2/issues/424


# TODO: change to dataframe
filter.log <- matrix(ncol = 2)

# TODO: check fastqPairedFilter
# https://github.com/benjjneb/dada2/issues/311

tic()
for(i in seq_along(fnFs)) {
  print(paste("Filering and Trimming sample: ", i))
  print(fnFs[[i]])
  print(fnRs[[i]])
  out <- dada2::filterAndTrim( 
                        fwd=fnFs[[i]], filt=filtFs[[i]],
                        rev=fnRs[[i]], filt.rev=filtRs[[i]],
                        trimLeft=trimLeft, trimRight=trimRight,
                        truncLen=truncLen,  # discard reads smaller then that and cut the rest
                        maxEE=maxEE, maxN=0, truncQ=QUALITY_THRESHOLD,  rm.phix=TRUE,
                        compress=FALSE, verbose=TRUE, multithread=TRUE
  )
  print(out)
  
  rownames(out) <- sample.names[i]  # change the name of the sample to a standard name like "34Sat2"
  filter.log <- rbind(filter.log, out)  # save filtered reads number
}
cat('Quality Filterin time:')
toc()  # 5839sec = 1.5 hours

# add a merger field to keep
filter.log <- cbind(filter.log, NA)
colnames(filter.log) <- c("reads.in", "reads.out", "merged")

########### Quality after filtering
report.path.filt <- paste0(filt_path.suf, "/fastQC")

# run all FASQC reports and save them to do Multiqc later
fastqcr::fastqc(fq.dir = filt_path.suf, # FASTQ files directory
                qc.dir = report.path.filt, # Results direcory
                threads = 5,                    # Number of threads
                fastqc.path = fastqc.path
)

# print the line for generation multi report
# multiqc /Users/alex/Projects_R/twins_microbiome_pipeline/data_set_bodyfl/fastq/no_primers/fastQC
print(paste("multiqc", report.path.filt))



########################################################################
###### BIG DATA sequential workflow in linear time #####################
if(length(filtFs) != length(filtRs)) stop("Forward and reverse files do not match!")


### LEARN ERROR RATES ###############
# https://github.com/benjjneb/dada2/issues/155
# for large data sets error rates should be estimated on subset of data - change to 40, here it is ok
tic()
errF <- learnErrors(filtFs, nbases=1e8, multithread = TRUE, randomize=TRUE, verbose=1, MAX_CONSIST=20)
errR <- learnErrors(filtRs, nbases=1e8, multithread = TRUE, randomize=TRUE, verbose=1, MAX_CONSIST=20)
cat("Error rate calculation time: ")
toc() # 9806.114 sec / 2.7h = now 2509sec

## plot error rates for control
plotErrors(errF)
plotErrors(errR)
save(errF, errR, file=file.path(files_intermediate_dada, dada.err.file)) 


### Big Data dada2 pipeline
# do it in a loop because we have a huge dataset (BigData workflow)
mergers <- list()
counter <- 0

tic()
for (sam in sample.names) {
  ### DEREPLICATION 
  # Dereplication is the process where all of the quality-filtered sequences are collapsed 
  # into a set of unique reads, which are then clustered into OTUs.
  derepF <- derepFastq(filtFs[[sam]])
  derepR <- derepFastq(filtRs[[sam]])
  
  ### SAMPLE INFERENCE 
  dadaF <- dada(derepF, err=errF, pool=TRUE, multithread = TRUE)
  dadaR <- dada(derepR, err=errR, pool=TRUE, multithread = TRUE)
  
  ### Merge forward and reverse reads
  # NOTE if merger is zero, it means your reads aren't overlapping after truncation - check trimming
  # mergePairs requires 20 nts of overlap by default
  # https://github.com/benjjneb/dada2/issues/419
  merger <- dada2::mergePairs(dadaF, derepF, dadaR, derepR)
  if (length(merger$sequence)==0){
    print(" !!!!!  You have a PROBLEM !!!!!! ")
    print("  Forward and Revers reads are not overlapping during merging! Please check trimming parameters!")
    print(" justConcatenate=TRUE has been used for this sample")
    merger <- dada2::mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE)  
    }
  mergers[[sam]] <- merger
  
  # add number of merged sequences to a filter.log
  filter.log[sam,"merged"] <- length(merger$sequence)
  
  counter <- counter + 1
  print("SAMPLE #", counter, "...", sam, ",  ", length(merger$sequence),  " merged sequences... Done.")
  print("----------")
}
cat("Total time of sample inference: ")
toc() # 8669sec=2.5h
save(mergers, file=file.path(files_intermediate_dada, mergers.file))


### SEQTAB: constructs a sequence table (analogous to an OTU table) from the list of samples.
seqtab.all <- dada2::makeSequenceTable(mergers, orderBy = "abundance")


### remove chimeras and save results
seqtab <- dada2::removeBimeraDenovo(seqtab.all, verbose = TRUE)
samples.names = rownames(seqtab)  # need it for futher Python data analysis

# named vector of all sequences
asv_sequences <- dada2::getSequences(seqtab)
prefix <- "asv_number"
suffix <- seq(1:length(asv_sequences))
asv.short.names <- paste(prefix, suffix, sep='_')  
names(asv_sequences) <- asv.short.names

save(seqtab, samples.names, asv_sequences, filter.log, file=file.path(files_intermediate_dada, seqtab.file)) 

print(" >>>  DADA2 ASV inference has been finished!")



# TODO: save seq tav as a CSV? OTU? any other format?
