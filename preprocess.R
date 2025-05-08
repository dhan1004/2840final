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
library(readr)
library(dbplyr)

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

sig_snps <- read_csv('significant_snps.csv',show_col_types = FALSE)

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
ld_matrix <- readRDS("ld_matrix.rds")

ld_df <- as.data.frame(as.table(ld_matrix))
colnames(ld_df) <- c("SNP1", "SNP2", "R2")
high_ld_pairs <- ld_df %>% filter(R2 > 0.8 & SNP1 != SNP2)
high_ld_pairs$SNP1 <- gsub("_.*$", "", high_ld_pairs$SNP1)
high_ld_pairs$SNP2 <- gsub("_.*$", "", high_ld_pairs$SNP2)
unique_snps <- unique(c(high_ld_pairs$SNP1, high_ld_pairs$SNP2))

# ensembl mapping?
ensembl = useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snp_gene_map <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start",
                 "chrom_end", "allele", "allele_1",
                 "minor_allele_freq", "ensembl_gene_stable_id", 
                 "ensembl_gene_name"),
  filters = "snp_filter",
  values = high_ld_pairs,
  mart = ensembl
)

# merge ld with gene mapping
ld_genes1 <- merge(
  high_ld_pairs,
  snp_gene_map,
  by.x = "SNP1",
  by.y = "refsnp_id",
  all.x = TRUE
)
ld_genes1_filtered <- ld_genes1[!is.na(ld_genes1$ensembl_gene_stable_id), ]

ld_genes2 <- merge(
  ld_genes1, 
  snp_gene_map, 
  by.x = "SNP2", 
  by.y = "refsnp_id", 
  all.x = TRUE)
colnames(ld_genes2) <- 
  c("SNP2", "SNP1", "R2", "CHR2", "CHR_START_2", "CHR_END_2", 
    "ALLELE_2_VAR", "ALLELE_2", "MAF_2", "GENE_STABLE_2", "GENE_NAME_2", "CHR1",
    "CHR_START_1", "CHR_END_1", 
    "ALLELE_1_VAR", "ALLELE_1", "MAF_1", "GENE_STABLE_1", "GENE_NAME_1")
ld_genes2_filtered <- 
  ld_genes2[!is.na(ld_genes2$GENE_STABLE_1) & !is.na(ld_genes2$GENE_STABLE_2) &
            ld_genes2$GENE_NAME_1 != "" & ld_genes2$GENE_NAME_2 != "", ]

gene_clusters <- ld_genes2_filtered %>%
  group_by(GENE_STABLE_1) %>%
  summarise(Connected_Genes = list(unique(GENE_STABLE_2)))

write.csv(ld_genes2_filtered, "mapped_ld_snps_genes.csv", row.names = TRUE)
saveRDS(ld_matrix, "mapped_ld_snps_genes.rds")
