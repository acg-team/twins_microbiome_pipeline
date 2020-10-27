# compare different taxonoby training sets on dada2 naive Bayes classifier

################# MSA, Tree and Taxonomy parameters
tools_param <- vector(mode="list", length=3)
names(tools_param) <- c("MSA_aligner", "tree_method", "tax_db", "tax_db_name")
 

##################################################################


tax.db.list <- c("silva/silva_nr99_v138_train_set.fa.gz",
                 "green_genes/gg_13_8_train_set_97.fa.gz", 
                 "rdp/rdp_train_set_16.fa.gz" 
                 )


# run dada2 naive Bayes classifier on all 3 databases
for(tax.db in seq_along(tax.db.list)) {
  tools_param$taxonomy_db <- tax.db
  tools_param$tax_db_name <- substring(tax.db, 1, 3)
  
  source("src/pipeline_dada2/5_Tax_Assign_dada2_RDP.R")
  
}