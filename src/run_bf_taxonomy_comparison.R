## compare taxonomy assignments done by DADA2 with three  taxomony DB: SILVA< GreenGenes, RDP


################# MSA, Tree and Taxonomy parameters
tools_param <- vector(mode="list", length=4)
names(tools_param) <- c("MSA_aligner", "tree_method", "tax_db", "tax_db_name")
 

##################################################################


tax.db.list <- c("silva/silva_nr99_v138_train_set.fa.gz",
                 "green_genes/gg_13_8_train_set_97.fa.gz", 
                 "rdp/rdp_train_set_16.fa.gz" 
                 )


##### run dada2 naive Bayes classifier on all 3 databases
for(tax.db in tax.db.list) {
  tools_param$tax_db <- tax.db
  tools_param$tax_db_name <- substring(tax.db, 1, 3)
  source("src/pipeline_dada2/5_Tax_Assign_dada2_RDP.R")
}

##### report all ASV which have different taxonomy assignments
# read all taxtables
load(file=file.path(files_intermediate_dada, "rdp_taxtab.RData")) 
tax.rdp <- taxtab

load(file=file.path(files_intermediate_dada, "sil_taxtab.RData")) 
tax.silva <- taxtab

load(file=file.path(files_intermediate_dada, "gre_taxtab.RData")) 
tax.green <- taxtab

# run comparison
mismatch.numbers <- 0
for(seq in rownames(tax.rdp)){
  # select taxonomy assignment for a given sequence and keep only Family/Genus
  rdp <- tax.rdp[seq,][c('Family','Genus')]
  silva <- tax.silva[seq,][c('Family','Genus')]
  green <- tax.green[seq,][c('Family','Genus')]
  
  rdp[is.na(rdp)] <- 0
  silva[is.na(silva)] <- 0
  green[is.na(green)] <- 0
  green <- replace(green, green=='', 0)
  
  # compare all three
  matches <- rdp == silva & rdp == green & silva == green

  
  if(!all(matches)) {
    mismatch.numbers <- mismatch.numbers + 1
    print(paste("a mismatch has been found for sequence: ", seq))
    print(rdp)
    print(silva)
    print(green)
  }
}

print(paste("Total mismatch number: ", mismatch.numbers))

# TODO: Create a csv file with report

