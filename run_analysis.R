# here you have to source all files like that
# NOTE: do not forget to set CONF variable with correct server/dataset values

# full workflow (1-2 days ona  a server)
source("src/analysis_twin/1_metadata.R")
source("src/2_file_names_parsing.R")

print("==================> long dada2 analysis has started...")
source("src/analysis_twin/3_BIG_dada_SV_table.R")

print("==================> Taxonomy assignment has started...")
source("src/analysis_twin/4_Tax_Assign.R")

print("==================> Phylogeny reconstraction has started...")
source("src/analysis_twin/5_Phylogeny.R")

source("src/analysis_twin/6_Create_Phyloseq_obj.R")


# Now PhyloSeq object is created and you can run analysis
