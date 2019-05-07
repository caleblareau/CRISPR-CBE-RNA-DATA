library(seqinr)

split_for_deeplift <- function(sample, n = 1000){
  editor <- "CBE"
  
  #Import sequence fasta files
  dir_seq <- paste0("../fastas/", editor, "-sequence")
  
  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)][grepl("chr21|chr22", seq_fastas)]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  seqs <- sapply(fasta_input, unname) %>% unlist() %>% unname()
  
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  
  # Split based on edit rate
  zero <- which(meta$editRate == 0 )
  mid <- which((meta$editRate >= 0.05) &  (meta$editRate < 0.20))
  high <- which(meta$editRate >= 0.5 )
  
  idx_zero <- zero[sample(1:length(zero), min(n, length(zero)))]
  idx_mid <- mid[sample(1:length(mid), min(n, length(mid)))]
  idx_high <- high[sample(1:length(high), min(n, length(high)))]
  
  # Export
  out_zero <- paste0("../fastas/CBE-forDeepLift/", sample, "-zero-sampled.fasta")
  out_mid <- paste0("../fastas/CBE-forDeepLift/", sample, "-mid-sampled.fasta")
  out_high <-   paste0("../fastas/CBE-forDeepLift/", sample, "-high-sampled.fasta")
  write.fasta(as.list(seqs[idx_zero]), seq_names[idx_zero],  out_zero, open = "w")
  write.fasta(as.list(seqs[idx_mid]), seq_names[idx_mid], out_mid, open = "w")
  write.fasta(as.list(seqs[idx_high]), seq_names[idx_high], out_high, open = "w")
  
  # Compress
  system(paste0("gzip ", out_zero))
  system(paste0("gzip ", out_mid))
  system(paste0("gzip ", out_high))
  
}
split_for_deeplift("89B")
split_for_deeplift("90B")
split_for_deeplift("160F")
split_for_deeplift("161F")
