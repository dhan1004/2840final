# install.packages("data.table")
# BiocManager::install(c("VariantAnnotation", "biomaRt", "GenomicRanges"))
# BiocManager::install("gwasvcf")
# remotes::install_github("mrcieu/gwasvcf")
# BiocManager::install("SeqArray")
# devtools::install_github("MRCIEU/ieugwasr")
# remotes::install_github("MRCIEU/genetics.binaRies")
# genetics.binaRies::get_plink_binary()
# install("BiocFileCache")

library(data.table)
library(VariantAnnotation)
library(biomaRt)
library(GenomicRanges)
library(gwasvcf)
library(SeqArray)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)
library(BiocManager)
library(rtracklayer)

# get significant SNPS ---------------------------------------------------------

vcf_file <- "C:/Users/sissy/coding/2840final/ieu-a-1239.vcf.gz"
vcf <- readVcf(vcf_file, "hg19")
vcf_index <- "C:/Users/sissy/coding/2840final/ieu-a-1239.vcf.gz.tbi"

lp_data <- assays(vcf)$LP
lp_data <- unlist(lp_data, use.names=FALSE)
pvals <- 10^(-(lp_data))

effect_size <- assays(vcf)$ES
effect_size <- unlist(effect_size, use.names=FALSE)

snp_ids <- rownames(vcf)
chr <- as.character(seqnames(vcf))
pos <- start(vcf)

# create data frame
if (length(snp_ids) == length(lp_data)) {
  # Create the data frame
  gwas_df <- data.frame(
    SNP = snp_ids,
    CHR = chr,
    POS = pos,
    LP = lp_data,
    P = pvals,
    ES = effect_size
  )
} else {
  stop("Mismatch in vector lengths. Check the data extraction process.")
}

sig_snps <- gwas_df %>% filter(P < 5e-8)
write.csv(sig_snps, "significant_snps.csv", row.names = FALSE)

# LD data ----------------------------------------------------------------------

# filter for top 10% of ES magnitude in sig snps
snp_list <- unique(sig_snps$SNP)
effect_size_threshold <- quantile(abs(sig_snps$ES), 0.9, na.rm = TRUE)
effect_size_threshold
top_snps <- sig_snps %>% filter(abs(ES) >= effect_size_threshold)
nrow(top_snps) # 3040 snps

snp_list <- as.list(top_snps$SNP)
str(snp_list)

plink_path <- genetics.binaRies::get_plink_binary()

# make ld matrix
ld_matrix <- ieugwasr::ld_matrix(
  top_snps$SNP,
  plink_bin = plink_path,
  bfile = "C:/Users/sissy/coding/2840final/reference_EUR/EUR"
)

write.csv(ld_matrix, "ld_matrix.csv", row.names = TRUE)
saveRDS(ld_matrix, "ld_matrix.rds")

ld_df <- as.data.frame(as.table(ld_matrix))
colnames(ld_df) <- c("SNP1", "SNP2", "R2")
high_ld_pairs <- ld_df %>% filter(R2 > 0.8 & SNP1 != SNP2)

# ucsc mapping? 
ld_snp_list <- unique(c(ld_df$SNP1, ld_df$SNP2))
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1000000, end = 1005000))
snp_track <- rtracklayer::ucscTableQuery("GRCh38.p14", range = gr)


# ensembl mapping?
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

