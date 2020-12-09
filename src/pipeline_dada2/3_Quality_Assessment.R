## Quality accessment of short reads 
# no primers


### LOAD PREVIOUS DATA
load(file=file.path(metadata_path, metadata.file)) 



### run QC QUALITY with fastqcr : ##############
library(fastqcr)

report.path <- paste0(raw_data_path, "/fastQC")

# run all reports
fastqcr::fastqc(fq.dir = raw_data_path, # FASTQ files directory
         qc.dir = report.path, # Results direcory
         threads = 5,                    # Number of threads
         fastqc.path = fastqc.path #"/Applications/BIOINF/FastQC.app/Contents/MacOS/fastqc"
)

# https://multiqc.info/
# you can run multiqc by  - 
# multiqc /Users/alex/Projects_R/twins_microbiome_pipeline/data_set_bodyfl/fastq/no_primers/fastQC

# aggregare the reports
qc <- qc_aggregate(report.path)

# https://cran.r-project.org/web/packages/fastqcr/readme/README.html
summary(qc)
qc_stats(qc)


