library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(reshape2)
library(stringr)
library(data.table)
library(diffloop)
library(seqinr)

"%ni%" <- Negate("%in%")

#---------------------------------------------------
# Essential functions for setting up all fasta files
#---------------------------------------------------

########
# Step 1
########

call_bam_readcount_ABE <- function(filepath_br, out_br){
  big_br_call <- paste0("zcat < ", filepath_br, " | awk '($3 == \"A\" || $3 == \"T\") && ($4 >=50) {print $0}'> ", out_br)
  system(big_br_call)
}


call_bam_readcount_CBE <- function(filepath_br, out_br){
  big_br_call <- paste0("zcat < ", filepath_br, " | awk '($3 == \"C\" || $3 == \"G\") && ($4 >=50) {print $0}'> ", out_br)
  system(big_br_call)
}

########
# Step 3
########
vienna_rna_bin <- "/data/aryee/caleb/deep/conda/envs/tf/bin/RNAfold"

secondaryStructureLaunch <- function(seq_fastas, library){
  
  # Set up shell script
  
  shell_file <- paste0("process_", library, "_structure.sh")
  
  
  shellbase <- c(
    "#!/bin/bash",
    "TEMP1", # #BSUB -J dl[1-3]
    "#BSUB -e /dev/null",
    "#BSUB -o /dev/null",
    "#BSUB -q normal",
    "",
    "# Pull input values",
    'TEMP2', 
    'input=$(cat $inputlist | awk -v ln=$LSB_JOBINDEX "NR==ln")',
    "",
    "# parse out individual variables from the string",
    "set -- $input",
    "old=$1",
    "new=$2",
    "",
    paste0(vienna_rna_bin, ' --noPS "${old}" | awk ', "'NR%3==0 || NR%3==1  {print $1}' > ", "${new}"),
    "gzip ${old}",
    "gzip ${new}",
    ""
  )
  
  # establish need
  need_structure_file <- paste0("need_structure/", library, "-structure-needArray.txt")
  needArrayLine <- paste0('inputlist="', need_structure_file, '"')
  
  # Define the structure file names
  new_structure_files <- gsub(".fasta$", ".structure.txt", gsub("-sequence", "-secondary", seq_fastas))
  
  # Establish the lookup / translation data frame
  need_df <- data.frame(
    seq_fastas, new_structure_files
  )
  need_number <- dim(need_df)[1]
  arrayHeaderLine <- paste0("#BSUB -J ",library,"[1-",as.character(need_number),"]")
  
  shellbase_here <- shellbase
  shellbase_here <- ifelse(shellbase_here == "TEMP1", arrayHeaderLine, shellbase_here)
  shellbase_here <- ifelse(shellbase_here == "TEMP2", needArrayLine, shellbase_here)
  
  # Export for bash array
  write.table(need_df, file = need_structure_file, 
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Export 
  write.table(data.frame(shellbase_here), file = shell_file, 
              sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Excute call
  #system(paste0("bsub < ", shell_file))
  shell_file
}



########
# Step 2
########

# Crude pass that pulls from gDNA -- does not consider any potential splicing
getSequenceSimple <- function(gr_query, reverse_complement = FALSE){
  
  gs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_query)
  
  # Flip the sequence if the user requests it
  if(reverse_complement){
    gs_rc <- reverseComplement(gs)
    gs <- gs_rc
  }
  
  return(as.character(gs))
}


getSequences_sensitive <- function(df, gtf, reverse_complement = FALSE, pad = 50){
  
  # Initialize a final results dataframe ; every should be straight forward but 'passed' is just essentially if
  # we could assemble a full sequence around the edit (unless we run out of exon real estate)
  df_out <- data.frame(
    chr = df$chr,
    start = df$start,
    total = df$TOTAL,
    strand = ifelse(reverse_complement, "R", "F"),
    gene = df$gene,
    editRate = round(df$editRate, 4),
    inOneTranscript = FALSE,
    passed = FALSE,
    sequence = "N", 
    originalIdx = 1:length(df$chr)
  )
  
  # Make a GRanges for both the padded and the SNV; keep only minimal sufficient meta data
  gr_pad <- padGRanges(makeGRangesFromDataFrame(df, keep.extra.columns = FALSE,
                                                seqnames.field = "chr", start.field = "start", end.field = "start"), pad = pad)
  gr_one <- makeGRangesFromDataFrame(df, keep.extra.columns = FALSE,
                                     seqnames.field = "chr", start.field = "start", end.field = "start")
  
  mcols(gr_one)$originalIdx <- df_out$originalIdx
  mcols(gr_one)$gene <- df_out$gene
  
  # See if the full padded sequence falls within the one transcript
  ov1 <- findOverlaps(gr_pad, gtf, type = "within")
  df_out$inOneTranscript <- 1:length(gr_pad) %in% queryHits(ov1)
  variants_gr_multi_onebp <- gr_one[1:length(gr_pad) %ni% queryHits(ov1)] # gr_one and gr_new have the same index order
  variants_gr_multi_pad <- gr_pad[1:length(gr_pad) %ni% queryHits(ov1)] # gr_one and gr_new have the same index order
  
  # Get the genomic sequences; we will have to update the multi-exon though
  df_out$sequence <- getSequenceSimple(gr_pad, reverse_complement = reverse_complement)
  
  # Remove instances where we observe an "N" in the reference genome (convenienty also when we couldn't define the sequence)
  df_out$passed <- (!grepl("N", df_out$sequence ) & nchar(df_out$sequence) == (2*pad + 1)) & df_out$inOneTranscript
  return(df_out)
}

# Function that 1) nominates variants and 2) annotates meta data for each library; these get written to a per-chromosome fasta
# for faster downstream computing

processFastaSample_CBE <- function(simple_name, step1_br_filter_file, dna_variants, gtf_forward, gtf_reverse, pad = 50){
  
  # Import data from bam readcount
  dt <- fread(step1_br_filter_file) %>% data.frame()
  colnames(dt) <- c("chr", "start", "REF", "TOTAL", "A", "C", "G", "T", "N")
  dt$poskey <- paste0(dt$chr, "_", dt$start)
  
  
  # Determine instances that overlap annotated DNA variants
  gr <- makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE,
                                 seqnames.field = "chr", start.field = "start", end.field = "start")
  ov_init <- findOverlaps(gr, dna_variants)
  
  # Find overlaps with GTF annotation + consider reference base
  pos_ov <- findOverlaps(gr, gtf_forward)
  neg_ov <- findOverlaps(gr, gtf_reverse)
  ref <- as.character(dt$REF)
  
  # Index over possibel variants to see ifi they meet critera
  ss_g <- 1:length(gr)
  possible_forward <- as.vector(ss_g %in% queryHits(pos_ov)  & ss_g %ni% queryHits(neg_ov) & ref == "C"  & seqnames(gr) %in% canonical)
  possible_reverse <- as.vector(ss_g %in% queryHits(neg_ov)  & ss_g %ni% queryHits(pos_ov) & ref == "G"  & seqnames(gr) %in% canonical)
  
  # Pull gene names // annotations for forward and reverse based on previous overlap
  gene_df_f <- data.frame(index = queryHits(pos_ov), gene = mcols(gtf_forward)$gene_name[subjectHits(pos_ov)]) %>% 
    group_by(index) %>% top_n(1, "gene") %>% data.frame()
  gene_forward_vec <- as.character(gene_df_f$gene); names(gene_forward_vec) <- as.character(gene_df_f$index)
  
  gene_df_n <- data.frame(index = queryHits(neg_ov), gene = mcols(gtf_reverse)$gene_name[subjectHits(neg_ov)]) %>% 
    group_by(index) %>% top_n(1, "gene") %>% data.frame()
  gene_reverse_vec <- as.character(gene_df_n$gene); names(gene_reverse_vec) <- as.character(gene_df_n$index)
  
  rm(gene_df_f); rm(gene_df_n); rm(pos_ov); rm(neg_ov)
  
  # Remove anything that was a reported DNA exome or common SNP
  possible_forward2 <- possible_forward & (1:length(gr) %ni% queryHits(ov_init))
  possible_reverse2 <- possible_reverse & (1:length(gr) %ni% queryHits(ov_init))
  
  # Subset and compute editing rates -- positive strand
  df_forward <- dt[possible_forward2,]
  df_forward$canonical_bases <- df_forward$`T` +  df_forward$C
  df_forward$editRate <- df_forward$`T`/( df_forward$TOTAL)
  df_forward$gene <- gene_forward_vec[as.character(which(possible_forward2))]
  
  # Subset and compute editing rates -- negative strand
  df_reverse <- dt[possible_reverse2,]
  df_reverse$canonical_bases <- df_reverse$A +  df_reverse$`G`
  df_reverse$editRate <- df_reverse$A/( df_reverse$TOTAL)
  df_reverse$gene <- gene_reverse_vec[as.character(which(possible_reverse2))]
  
  #rm(dt); rm(dna_variants)
  
  # A final filter to remove 1) bases that have lots of non-canonical bases and 2) that lack coverage at canonical bases
  df_forward %>% filter(canonical_bases >= 50 & (canonical_bases/TOTAL) >= 0.98) -> df_forward
  df_reverse %>% filter(canonical_bases >= 50 & (canonical_bases/TOTAL) >= 0.98) -> df_reverse
  
  # Define sequences and essential attributes per-edit
  df_forward_final <- getSequences_sensitive(df_forward, gtf_forward)
  df_reverse_final <- getSequences_sensitive(df_reverse, gtf_reverse, reverse_complement = TRUE)
  df_final <- rbind(df_forward_final, df_reverse_final)
  
  # Final annoying overlap to get annotation
  gr_last <- makeGRangesFromDataFrame(df_final, keep.extra.columns = FALSE,
                                      seqnames.field = "chr", start.field = "start", end.field = "start", strand.field = "")
  
  # Import UTR annotation
  x3utr <- readRDS("/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/hg38_annotation/x3utr.hg38.gr.rds"); ov3 <- findOverlaps(gr_last, x3utr)
  x5utr <- readRDS("/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/hg38_annotation/x5utr.hg38.gr.rds"); ov5 <- findOverlaps(gr_last, x5utr)
  xCDS <- readRDS("/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/hg38_annotation/cds.hg38.gr.rds"); ovcds <- findOverlaps(gr_last, xCDS)
  
  df_final$annotation <- ifelse(1:dim(df_final)[1] %in% queryHits(ov5),"5UTR",
                                ifelse(1:dim(df_final)[1] %in% queryHits(ov3), "3UTR",
                                       ifelse(1:dim(df_final)[1] %in% queryHits(ovcds), "CDS", "intron")))
  
  # Combine and export per chromosome (to faciliate easier i/o, etc.)
  possible_chrs <- as.character(unique(df_final$chr))
  
  lout <- lapply(possible_chrs, function(chr){
    
    # Get attributes
    df_chr <- df_final[df_final$chr == chr, ]
    positions <- paste0(as.character(df_chr$chr), "-", as.character(df_chr$start))
    strand <- as.character(df_chr$strand)
    rates <- as.character(df_chr$editRate)
    coverage <- as.character(df_chr$total)
    annotation <-  as.character(df_chr$annotation)
    gene <-  as.character(df_chr$gene)
    sequences <- as.character(df_chr$sequence)
    
    # Write to fasta
    fasta_file <- paste0("../fastas/CBE-sequence/", simple_name, "." ,chr, ".fasta")
    fasta_names <- paste0(simple_name, "_", positions, "_", strand, "_", rates, "_", coverage, "_",  annotation, "_", gene)
    write.fasta(as.list(sequences), fasta_names, fasta_file, open = "w")
    fasta_file
  })
  
  # Figure out where we lost potential variants as a quality control measure
  allVariants0 <- length(gr) # everything with a reference nucleotide that could be edited (either strand) and a raw 50X coverage
  specificNucleotides1 <- sum(possible_forward) + sum(possible_reverse) # positions where the reference nucleotide matches the strand for a specific gene (non-overlapping) AND is not a 
  nonDNAvariants2 <- sum(possible_forward2) + sum(possible_reverse) # Positions that did not overlap a known variant
  canonicalVariants3 <- dim(df_forward)[1] + dim(df_reverse)[1] 
  finalVariants4 <- sum(df_forward_final$passed) + sum(df_reverse_final$passed)
  
  # Write out to qc data frame
  qc_df <- data.frame(
    simple_name, 
    allVariants0,
    specificNucleotides1,
    nonDNAvariants2,
    canonicalVariants3,
    finalVariants4,
    forward_strand = sum(df_forward_final$passed),
    reverse_strand = sum(df_reverse_final$passed)
  )
  
  write.table(qc_df, file = paste0("../fastas/CBE-qc/", simple_name, "-qc.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  return(lout)
}
