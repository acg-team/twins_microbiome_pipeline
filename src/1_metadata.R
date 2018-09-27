# @Alex
# LOAD metadate of all samples from
#http://www.ebi.ac.uk/ena/data/view/PRJEB13747
# and put it into a file

install.packages("tictoc")
require(XML)
require(tictoc)

# create a path variable to access data and processed-filtered data
project_path <- "~/Projects_R/twins_microbiome_pipeline"
setwd(project_path)
source("src/load_initialize.R")

#####: Extract metadata desription information from different files
metadata.raw <- read.table(file.path(data_path,"PRJEB13747.txt"))

# initialize vectors of different features
age_vector     <- c()
twin_id_vector <- c()
collection_date_vector <- c()
sex_vector       <- c()
zygosity_vector  <- c()
family_id_vector <- c()
host.subject.id_vector <- c()

sample_number = 3289 ##3289

# extract metainformation for each sample from anXML at ebi
for (sample_idx in 2:sample_number){ 
  # construct url
  # http://www.ebi.ac.uk/ena/data/view/ERS1131064&display=xml&download=xml&filename=ERS1131064.xml
  xml.url <- paste("http://www.ebi.ac.uk/ena/data/view/", metadata.raw[sample_idx,"V3"], "&display=xml&download=xml&filename=",metadata.raw[sample_idx,"V3"],".xml", sep="" )
  xmldata <- xmlParse(xml.url)
  
  collection_date <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="collection_date"]/VALUE')[[1]] )
  age       <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host age"]/VALUE')[[1]] )
  twin_id   <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="twin_id"]/VALUE')[[1]] )
  sex       <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host sex"]/VALUE')[[1]] )
  zygosity  <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="zygosity"]/VALUE')[[1]] )
  family_id <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="family_id"]/VALUE')[[1]] )
  host.subject.id <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host subject id"]/VALUE')[[1]] )
  
  collection_date_vector <- c(collection_date_vector,collection_date )
  age_vector       <- c(age_vector, age)
  twin_id_vector   <- c(twin_id_vector, twin_id)
  sex_vector       <- c(sex_vector, sex)
  zygosity_vector  <- c(zygosity_vector, zygosity)
  family_id_vector <- c(family_id_vector, family_id)
  host.subject.id_vector <- c(host.subject.id_vector, host.subject.id)
  
}

df.metadata <- data.frame(
  file            = metadata.raw[2:sample_number,"V4"],
  collection_date = collection_date_vector,
  twin_id         = twin_id_vector,
  sex             = sex_vector,
  zygosity        = zygosity_vector,
  family_id       = family_id_vector,
  host.subject.id = host.subject.id_vector,
  age             = age_vector
)

df.metadata.ordered <- df.metadata[ order(df.metadata["family_id"]), ]

# cache to disk to avoid downloading again
save(df.metadata, df.metadata.ordered, file=file.path(result_path, "metadata.RData")) 

# export to human readable format
write.table(df.metadata, file.path(result_path,"metadata.csv"), sep=",")
write.table(df.metadata.ordered, file.path(result_path,"metadata_ordered.csv"), sep=",")


