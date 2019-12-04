# UK twins microbiome processing pipeline
Bioinformatics workflow for processing raw Illumina reads of microbiome.

INPUTS: raw illimina data, metainformation (age, sex etc)
OUTPUT: exact sequence variant table (SV table), taxomony assignments, phylogeny of samples

Contains a complete workflow from extracting metainformation to assigning taxa information and phylogenetic tree of beta diversion.
Pipeline as from https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html


## Raw Data
- download all raw illumina reads from http://www.ebi.ac.uk/ena/data/view/PRJEB13747
- place them into /data/raw

## WorkFlow: run subsequently
### 1. [ 1_metadata.R  ]

A result od this script should be creating a dataframe with metadata [ twin_id/ sex/ zigosity/ etc ] attached to each sample's name

  - downloads meta information attatched to every sample, age, family_id etc from http://www.ebi.ac.uk/ena/data/view/ERS1131064&display=xml&download=xml&filename=ERS1131064.xml
  - parse XML into dataframe
  - save in into data/metadata/metadata.RData.
Possible inprovements: If nessesary, modify the script to get more features.

### 2. [ 2_BIG_dada_SV_table.R ]

WARNING: most computationally demanding script! ~1 day on 4 core server with 32G memory

This script inplements dada2 pipeline to convert RAW Illumina reads into Single Variants (SV) table (analogous to OTU, but on variant level, not species).

The pipeline is desribed here - https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

RESULTS: the result is a seqtab abandance matrix (3288 samples[ERR138, ...] x 8299 sequence variants) with abundancies in every cell.
saved into "seqtab_q15.RData" (q means trim parameter used during run)

IMPROMEMENTS:
  - play with trim parameters for raw reads
  - play wit other quality prameyers
  
  
### 3. [ 3_Tax_Assign.R ]

For each sequence variant deduced during previous step (8299) assign a taxomomy
RESULT: taxtab_g15.RData file, [ 8299 sequence variants x [Kingdom, Phulym, etc]]

IMPROVEMENTS:
   - silva assignment is not very accurate, check another 
   
   
### 4. [ 4_Phylogeny.R ]

Create a phylogeny out of all deduced sequence 8299 variants as follows
  - MSA with Muscle or ClustalW (seems very slow, need to investigate its advantages if any)
  - create a guide tree with NJ
  - build a tree with Phangorn, and parameters : model="GTR", optInv=TRUE, optGamma=TRUE

WARNING: has to be improved

POSSIBLE IMPROVEMENTS:
  - try PhyML ;
  - since we have short sequences (250) we might try to use phylogenies with indel in Bayesian settings;
  - read about metagenomic taxonomy and functional assignment;
  
  
### 5. [ 5_Create_Phyloseq_obj.R ]

Create and save Phyloseq object for further manipulation and visulization of microbiome data
https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#aim


### 6. [ 5_Phyloseq_Analysis.R ]
Use Phyloseq to 
  - calculate a distance matrix btw all 3288 samples based on abundance/ Tree information (UNIFRAC)
  - plot hierarsicul clustering
  - try to do PCA / PCoA to detect variations



## Possible futher directions
#### Calculate the metric based on taxonomy  
  - Phylogenetics: use PhyML, then try reconstruction with indels and TKF model
  - generate UNIFRAC metrics (needs a good phylo-tree first)
  
#### Calculate the metric based on abandance
  - do we already have such a metric?
  - can we combine UNIFRAC with abandancies?
  
#### Other approaches  
  - check vusualization with Phyloseq
  - figure out the best parameters for quality filtering in DADA2 pipeline
  - generate OTU as well as SV table
  
  
## References and packages
#### packages to study
- ape - Analyses of Phylogenetics and Evolution
- RAxML - create a phylo tree for big abount of short sequnces (16S)
- phangorn - phylo tree building package
- Phyloseq - package to analyse microbial communities


#### Articles
- 


