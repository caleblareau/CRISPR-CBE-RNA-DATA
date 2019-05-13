
# Classify variants
dir <- "/data/joung/caleb/base_editing/sequence-determinants/sequence-RNA-editing/data/hg38_annotation/"
x5utr <- readRDS(paste0(dir, "x5utr.hg38.gr.rds"))
x3utr <- readRDS(paste0(dir, "x3utr.hg38.gr.rds"))
exons <- readRDS(paste0(dir, "exons.hg38.gr.rds"))
introns <- readRDS(paste0(dir, "introns.hg38.gr.rds"))

annotate_variants <- function(vars_gr){
  ov_x5utr <- queryHits(findOverlaps(vars_gr, x5utr))
  ov_x3utr <- queryHits(findOverlaps(vars_gr, x3utr))
  ov_exons <- queryHits(findOverlaps(vars_gr, exons))
  ov_introns <- queryHits(findOverlaps(vars_gr, introns))
  pos <- 1:length(vars_gr)
  
  anno <- ifelse(pos %in% ov_x5utr, "5UTR", 
                 ifelse(pos %in% ov_x3utr, "3UTR", 
                        ifelse(pos %in% ov_exons, "exon", 
                               ifelse(pos %in% ov_introns, "intron", "unknown"))))
}


# Function to yield True / False vector of length(anno2) such that the proprtions match those in anno1
optimize_matched_subset <- function(anno1, anno2){
  
  stopifnot(length(anno2) > length(anno1))
  
  props <- table(anno1)/length(anno1) %>% as.vector()
  
  # Figure out which is under represented, proprtionally, in anno2
  expectation <- length(anno2) *props
  observed <- table(anno2)
  ratio <- expectation/observed
  lagging <- which(max(ratio) == ratio) %>% names()
  
  # Determine number of observations that can go in each group and be balanced
  n_total <- round(observed[lagging]/props[lagging]) %>% unname() - 10 # total number that we can have // safety
  new_target <- round(n_total * props)
  
  # Verify that that we can fill each group
  stopifnot(all(sapply(1:length(new_target), function(i) observed[i] >= new_target[i])))
  
  # Do a subsample for each class
  classes <- names(new_target)
  index <- 1:length(anno2)
  set.seed(100)
  sapply(classes, function(class){
    class_options <- index[anno2 == class]
    sample(class_options,  size = new_target[class])
  }) %>% unlist() %>% unname() %>% sort() -> keepers
  return(index %in% keepers)
}
