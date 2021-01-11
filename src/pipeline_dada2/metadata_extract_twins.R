# @AlexY, 2020
# NOTE: only for twin dataset!
# LOAD metadata of all samples from EBI
# http://www.ebi.ac.uk/ena/data/view/PRJEB13747
# and save it into a CSV file

# it is only neccesary to run it once, then files will be cached on disk
# NOTE: loar.R and configure.R must be already loaded in the outside pipeline.R file!
library(httr)


#####: Extract metadata desription information from different files
metadata.raw <- read.table(file.path(metadata_path,"PRJEB13747.txt"))

# initialize vectors of different features
age_vector     <- c()
twin_id_vector <- c()
collection_date_vector <- c()
sex_vector       <- c()
zygosity_vector  <- c()
family_id_vector <- c()
host.subject.id_vector <- c()

sample_number = 3289 ##3289

# extract meta information from anXML / ebi for each sample
for (sample_idx in 2:sample_number){ 
  # construct url
  # http://www.ebi.ac.uk/ena/data/view/ERS1131064&display=xml&download=xml&filename=ERS1131064.xml
  xml.url <- paste("http://www.ebi.ac.uk/ena/data/view/", metadata.raw[sample_idx,"V3"], "&display=xml&download=xml&filename=",metadata.raw[sample_idx,"V3"],".xml", sep="" )
  #xmldata <- xmlParse(xml.url)
  xmldata <- xmlParse(rawToChar(GET(xml.url)$content))
  
  collection_date <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="collection_date"]/VALUE')[[1]] )
  age       <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host age"]/VALUE')[[1]] )
  twin_id   <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="twin_id"]/VALUE')[[1]] )
  sex       <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host sex"]/VALUE')[[1]] )
  zygosity  <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="zygosity"]/VALUE')[[1]] )
  family_id <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="family_id"]/VALUE')[[1]] )
  host.subject.id <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="host subject id"]/VALUE')[[1]] )
  country <- xmlValue(  getNodeSet(xmldata,'//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG="geographic location (country and/or sea)"]/VALUE')[[1]] )
  
  collection_date_vector <- c(collection_date_vector,collection_date )
  age_vector       <- c(age_vector, age)
  twin_id_vector   <- c(twin_id_vector, twin_id)
  sex_vector       <- c(sex_vector, sex)
  zygosity_vector  <- c(zygosity_vector, zygosity)
  family_id_vector <- c(family_id_vector, family_id)
  host.subject.id_vector <- c(host.subject.id_vector, host.subject.id)
  country_vector <- c(country_vector, country)
}

df.metadata <- data.frame(
  file            = metadata.raw[2:sample_number,"V4"],
  collection_date = collection_date_vector,
  twin_id         = twin_id_vector,
  sex             = sex_vector,
  zygosity        = zygosity_vector,
  family_id       = family_id_vector,
  host.subject.id = host.subject.id_vector,
  age             = age_vector,
  country         = country_vector
)
rownames(df.metadata) <- df.metadata$file # use ERR as row names instead of 1,2,3...

df.metadata.ordered <- df.metadata[ order(df.metadata["family_id"]), ]  # order by family_id

#### filter the list of samples to have 2 samples for eash twin
# get family list with >3 members
number_of_samples_in_family <- table(df.metadata.ordered$family_id)
family_id_of_large_families <- rownames(number_of_samples_in_family[number_of_samples_in_family>3])

# get twin list with more then 1 samples
number_of_samples_of_twin <- table(df.metadata.ordered$twin_id)
twin_id_more_1sample <- rownames(number_of_samples_of_twin[number_of_samples_of_twin>1])

#filter out twins with less then 2 samples
df.metadata.4timepoints <- df.metadata.ordered[df.metadata.ordered$twin_id %in% twin_id_more_1sample,]

# filter out families with lt 4 samples
df.metadata.4timepoints <- df.metadata.4timepoints[df.metadata.4timepoints$family_id %in% family_id_of_large_families,]


# cache to disk to avoid downloading again
save(df.metadata, df.metadata.ordered, df.metadata.4timepoints, file=file.path(metadata_path, metadata.file)) 

# export to human readable format
write.table(df.metadata, file.path(metadata_path,"metadata.csv"), sep=",")
write.table(df.metadata.ordered, file.path(metadata_path,"metadata_ordered.csv"), sep=",")


