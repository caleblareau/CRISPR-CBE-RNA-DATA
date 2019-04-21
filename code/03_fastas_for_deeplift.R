library(seqinr)

split_for_deeplift <- function(sample, n = 1000, chrs = 1:2){
  editor <- "CBE"
  # Import sequence fasta files
  dir_seq <- paste0("../fastas/", editor, "-sequence")

  seq_fastas <- list.files(dir_seq, pattern = paste0(sample), full.names = TRUE)
  seq_fastas <- seq_fastas[grepl(".gz", seq_fastas)][chrs]
  
  # Import data into a list
  fasta_input <- sapply(seq_fastas, function(faf) unlist(read.fasta(faf, as.string = TRUE)))
  seq_names <- sapply(fasta_input, names) %>% unlist() %>% unname()
  seqs <- sapply(fasta_input, unname) %>% unlist() %>% unname()
  
  meta <- data.frame(stringr::str_split_fixed(seq_names, pattern = "_", n = 7), stringsAsFactors = FALSE)
  colnames(meta) <- c("Library", "chr_pos", "strand", "editRate", "coverage",  "annotation", "gene")
  
  # Split based on edit rate
  zero <- which(meta$editRate == 0 )
  mid <- which(meta$editRate >= 0.05 &  meta$editRate < 0.20)
  high <- which(meta$editRate >= 0.5 )
  
  idx_zero <- zero[sample(1:length(zero), min(n, length(zero)))]
  idx_mid <- mid[sample(1:length(mid), min(n, length(mid)))]
  idx_high <- high[sample(1:length(high), min(n, length(high)))]
  
  # Export
  write.fasta(as.list(seqs[idx_zero]), seq_names[idx_zero],
              paste0("../fastas/CBE-forDeepLift/", sample, "-zero-sampled.fasta"), open = "w")
  
  write.fasta(as.list(seqs[idx_mid]), seq_names[idx_mid],
              paste0("../fastas/CBE-forDeepLift/", sample, "-mid-sampled.fasta"), open = "w")
  
  write.fasta(as.list(seqs[idx_high]), seq_names[idx_high],
              paste0("../fastas/CBE-forDeepLift/", sample, "-high-sampled.fasta"), open = "w")

}
split_for_deeplift("89B")
