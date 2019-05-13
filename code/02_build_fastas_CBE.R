library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(reshape2)
library(stringr)
library(data.table)
library(diffloop)
library(seqinr)

# Import annotation data only if it is missing from the environment
if(!exists("exome_hepg2")){
  source("00_helperdata.R")
}

# Import functions to facilitate data preparation
source("01_essentialFunctions-cbesecure.R")

do_all_CBE <- function(simple_name, full_bamreadcount_library_file, cell_type = "HEK293"){
  
  if(cell_type == "HEK293"){
    dna_variants=exome_hek293t
  }else if(cell_type == "HepG2"){
    dna_variants=exome_hepg2
  }
  
  # Create new variables based on input
  step1_br_filter_file <- paste0("bam-readcount/" , simple_name, "-HQcounts.tsv")
  
  # Step 1 - Determine universe of possible edits / non-edits from bam-readcount data
  call_bam_readcount_CBE(full_bamreadcount_library_file , step1_br_filter_file)
  
  # Step 2 - Determine universe of possible edits / non-edits from bam-readcount data
  new_fastas <- processFastaSample_CBE(simple_name, step1_br_filter_file, dna_variants, gtf_forward, gtf_reverse, pad = 50)
  
  # Step 3 - Annotate secondary structure for all files generated in step 2
  new_structure_files <- secondaryStructureLaunch(unlist(new_fastas), simple_name)
  
  # Step 4 - remove/compress for economy
  system(paste0("rm ", step1_br_filter_file))
  simple_name
}


if(TRUE){
  
  #do_all_CBE(simple_name = "162B",
  #           full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/162B.all_positions.txt.gz",
  #           cell_type = "HEK293")
  
  do_all_CBE(simple_name = "162F",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/162F.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "162G",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/162G.all_positions.txt.gz",
             cell_type = "HEK293")
  
  do_all_CBE(simple_name = "162H",
             full_bamreadcount_library_file = "/data/joung/caleb/base_editing/exp7-hiseq/bam-readcount/162H.all_positions.txt.gz",
             cell_type = "HEK293")
  
}


