library(data.table)
library(diffloop)
library(seqinr)
library(rtracklayer)
library(MutationalPatterns)

# Import GTF and split by strand
gtf <- import("/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/gencode.v29.appris_principle_1.gtf")

# Filter for canonical chromosomes 
canonical <- paste0("chr", c(as.character(1:22), "X"))
gtf <- gtf[(seqnames(gtf) %in% canonical)]

# Split based on strand
gtf_forward <- gtf[strand(gtf) == "+"]
gtf_reverse<- gtf[strand(gtf) == "-"]
rm(gtf)

exome_hek293t_file <- "/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/dna_variants/HEK293_allDNA.rds"
exome_hepg2_file <- "/data/joung/caleb/base_editing/sequence-determinants/CRISPR-BE-RNA-editing/data/dna_variants/HepG2_allDNA.rds"

# Cache the variants needed for filtering
if(file.exists(exome_hek293t_file) &&file.exists(exome_hepg2_file)) {
  
  exome_hek293t <- readRDS(exome_hek293t_file)
  exome_hepg2 <- readRDS(exome_hepg2_file)
  
} else {
  
  # Import exome vcf ++ common variants -- not particularly elegant but fast and relatively conservative
  SNPs <- fread("/data/joung/caleb/base_editing/hg38/dbsnp_138.hg38.bed2") %>% data.frame() %>%
    makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V2")
  
  exome_hek293t <- c(fread("/data/joung/caleb/base_editing/exome/variantCalls/gatk3.8-erisone/allHEK293A.vcf", skip = 200, select = c("V1", "V2"), header = FALSE) %>%
                       data.frame() %>% filter(V1 %in% canonical) %>% data.frame() %>%
                       makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V2"), SNPs) %>% unique()
  
  exome_hepg2 <- c(fread("/data/joung/caleb/base_editing/exome/variantCalls/gatk3.8-erisone/allHepG2A.vcf" , skip = 200, select = c("V1", "V2"), header = FALSE) %>%
                     data.frame() %>% filter(V1 %in% canonical) %>% data.frame() %>%
                     makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V2"), SNPs) %>% unique()
  
  saveRDS(exome_hek293t, file = exome_hek293t_file)
  saveRDS(exome_hepg2, file = exome_hepg2_file)
  
}
