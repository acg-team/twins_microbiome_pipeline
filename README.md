# twins_microbiome_pipeline
Bioinformatics procession of raw Illumina reads for microbiome analysis of UK twins study with dada2 pipeline

Contains a complete workflow from extracting metainformation to assigning taxa information and phylogenetic tree of beta diversion.
Pipeline as from https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html


## raw Data and third-party data
- download all raw illumina reads from http://www.ebi.ac.uk/ena/data/view/PRJEB13747
- place them into /data/raw

## run subsequently
1. [ 1_metadata.R  ]

A result od this script should be creating a dataframe with metadata [ twin_id/ sex/ zigosity/ etc ] attached to each sample's name

  - downloads meta information attatched to every sample, age, family_id etc from http://www.ebi.ac.uk/ena/data/view/ERS1131064&display=xml&download=xml&filename=ERS1131064.xml
  - parse XML into dataframe
  - save in into data/metadata/metadata.RData.
Possible inprovements: If nessesary, modify the script to get more features.

2. [ 2_BIG_dada_SV_table.R ]

WARNING: most computationally demanding script! ~1 day on 4 core server with 32G memory

This script inplements dada2 pipeline to convert RAW Illumina reads into Single Variants (SV) table (analogous to OTU, but on variant level, not species).

The pipeline is desribed here - https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

RESULTS: the result is a seqtab abandance matrix (3288 samples[ERR138, ...] x 8299 sequence variants) with abundancies in every cell.
saved into "seqtab_q15.RData" (q means trim parameter used during run)

IMPROMEMENTS:
  - play with trim parameters for raw reads
  - play wit other quality prameyers
  
  
3. [ 3_Tax_Assign.R ]

For each sequence variant deduced during previous step (8299) assign a taxomomy
RESULT: taxtab_g15.RData file, [ 8299 sequence variants x [Kingdom, Phulym, etc]]

IMPROVEMENTS:
   - silva assignment is not very accurate, check another 
   
   
4. [ 4_Phylogeny.R ]

Create a phylogeny out of all deduced sequence 8299 variants as follows
  - MSA with ClustalW
  - create a guide tree with NJ
  - build a tree with Phangorn, and parameters : model="GTR", optInv=TRUE, optGamma=TRUE

WARNING: has to be improved

POSSIBLE IMPROVEMENTS:
  - try PhyML ;
  - since we have short sequences (250) we might try to use phylogenies with indel in Bayesian settings;
  - read about metagenomic taxonomy and functional assignment;
  
  
5. [ 5_Create_Phyloseq_obj.R ]

Create  Phyloseq object for further manipulation and vusulization of microbiome data
https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#aim

TODO:
  - study Phyloseq


