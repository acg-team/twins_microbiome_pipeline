# UK twins microbiome processing pipeline
Bioinformatics pipeline for processing 16S amplican raw Illumina reads and subsequent microbiome analysis.

INPUTS: raw illimina data, metainformation (age, sex etc)
OUTPUT: 
- exact sequence variant table (SV table - dada2)
- taxomony assignments (silva, dada2)
- phylogeny of samples (RAxML)
- combine all in Phyloseq object
- re-done everythong with QIIME2 for OTU

Contains a complete workflow from extracting metainformation to assigning taxa information and phylogenetic tree of beta diversion.
Pipeline as from https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html


## Raw Data
- download all raw illumina reads from http://www.ebi.ac.uk/ena/data/view/PRJEB13747
- place them into /data_set_twin/raw

## Pepiline: run subsequently
### 0. [load.R  / configure.R]
Those scripts load R packages and previously processed data. No need to run them separatelly, they are included in each script below.
In future, that structure might be substituted with NextFlow framework


### 1. [ src/pipeline_dada2/1_metadata.R  ]

A result od this script should be creating a dataframe with metadata [ twin_id/ sex/ zigosity/ etc ] attached to each sample's name

  - downloads meta information attatched to every sample, age, family_id etc from http://www.ebi.ac.uk/ena/data/view/ERS1131064&display=xml&download=xml&filename=ERS1131064.xml
  - parse XML into dataframe
  - save in into data/metadata/metadata.RData.
Possible inprovements: If nessesary, modify the script to get more features.

### 2. [ src/pipeline_dada2/3_BIG_dada_SV_table.R ]

WARNING: most computationally demanding script! ~1 day on 4 core server with 32G memory

This script inplements dada2 pipeline to convert RAW Illumina reads into Single Variants (SV) table (analogous to OTU, but on variant level, not species).

The pipeline is desribed here - https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

RESULTS: the result is a seqtab abandance matrix (3288 samples[ERR138, ...] x 8299 sequence variants) with abundancies in every cell.
saved into "seqtab_q15.RData" (q means trim parameter used during run)

IMPROMEMENTS:
  - play with trim parameters for raw reads
  - play with other quality prameyers
  
  
### 3. [ src/pipeline_dada2/4_Tax_Assign.R ]

For each sequence variant deduced during previous step (8299) assign a taxomomy
RESULT: taxtab_g15.RData file, [ 8299 sequence variants x [Kingdom, Phulym, etc]]

IMPROVEMENTS:
   - silva assignment is not very accurate, check another 
   
   
### 4. [ src/pipeline_dada2/5_Phylogeny.R ]

Create a phylogeny out of all deduced sequence 8299 variants as follows
  - MSA with Muscle or ClustalW (seems very slow, need to investigate its advantages if any)
  - create a guide tree with NJ
  - build a tree with Phangorn, and parameters : model="GTR", optInv=TRUE, optGamma=TRUE (fast option)
  - Run  RAxML (slow, but accurate option) - need to be installed on local machine!
  
  
### 5. [ src/pipeline_dada2/6_Create_Phyloseq_obj.R ]

Combine everithing in one place.
Create and save Phyloseq object for further manipulation and visulization of microbiome data
https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#aim


### [ src/analysis_twin/ ]




## Possible futher directions

  
  
## References and packages
#### packages to study
- ape - Analyses of Phylogenetics and Evolution
- RAxML - create a phylo tree for big abount of short sequnces (16S)
- phangorn - phylo tree building package
- Phyloseq - package to analyse microbial communities


#### Articles
- 


