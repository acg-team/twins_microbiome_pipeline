# 1. read the RAW data folder
# 2. filter them if neccesary
# 3. split then into two lists with forward and reverse reads
# by the and of this script we have fnFs / fnRs as lists of reads names


load(file=file.path(metadata_path, metadata.file)) 


###### create lists of file names of forward reads (fnFs) and reverse reads (fnRs) ######

# get all sample's file names to be processed by scanning the raw data folder
fns <- sort(list.files(raw_data_path, full.names = TRUE))
print(paste("Samples (reverse and forward) left BEFORE filtering::", length(fns)))


if(conf$dataset == 'TWIN'){
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
  
  save(df.metadata, df.metadata.ordered, df.metadata.4timepoints, sample.names, fnFs, fnRs, file=file.path(metadata_path, metadata.file)) 

} else if(conf$dataset == 'BFL'){
  # split all filenames into Forward and Reverse reads files
  fnFs <- fns[grepl(".R1.fastq", fns)]
  fnRs <- fns[grepl(".R2.fastq", fns)]
  
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
  
  save(df.metadata, sample.names, fnFs, fnRs, file=file.path(metadata_path, metadata.file)) 
  
} else {
  stop("WRONG dataset configuration name!")
}


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


# Result:
head(fnFs)
head(fnRs)
head(sample.names)
