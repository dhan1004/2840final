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

if (length(snp_ids) == length(lp_data)) {
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

#Luke Mapping
mapped_ld_snp_genes <- readRDS("../../2840_data/mapped_ld_snps_genes.rds")
mapped_ld_snp_genes_csv <- read.csv("../../2840_data/mapped_ld_snps_genes.csv")
ld_matrix_csv <- read.csv("../../2840_data/ld_matrix.csv")
significant_snps_csv <- read.csv("../../2840_data/significant_snps.csv")

sig_snps_high_ES <- significant_snps_csv %>%
  dplyr::filter(
      abs(ES) >= 0.02
  )
sig_snps_ES <- sig_snps_high_ES$SNP
mapped_ld_snp_genes_csv_high_ES <- mapped_ld_snp_genes_csv %>% dplyr::filter(
  SNP1 %in% sig_snps_ES & SNP2 %in% sig_snps_ES)

#mapped_ld_snp_genes_csv_high_ES_high_R_unique_genes <- mapped_ld_snp_genes_csv_high_ES %>%
#  dplyr::filter(R2 > 0.8 & GENE_NAME_2 != GENE_NAME_1)


mapped_ld_snp_genes_csv_FINAL <- mapped_ld_snp_genes_csv_high_ES %>%
  dplyr::mutate(
    CHR_START_1 = as.numeric(CHR_START_1),
    CHR_START_2 = as.numeric(CHR_START_2)
  ) %>%
  dplyr::filter(
    !(CHR1 == CHR2 & abs(CHR_START_1 - CHR_START_2) <= 1e6)
  ) %>%
  dplyr::filter(
    R2 > 0.8
  )

sig_snps_ES <- sig_snps_high_ES$SNP
all_linked_SNPs <- unique(c(mapped_ld_snp_genes_csv_FINAL$SNP1,mapped_ld_snp_genes_csv_FINAL$SNP2))
all_linked_SNPs_DATA <- sig_snps_high_ES %>% dplyr::filter(SNP %in% all_linked_SNPs)

all_UNLINKED_SNPs_DATA <- sig_snps_high_ES %>% dplyr::filter(!SNP %in% all_linked_SNPs)

SNP_to_ensembl_second <- mapped_ld_snp_genes_csv %>%
  dplyr::select(SNP = SNP2, GENE_NAME = GENE_STABLE_2)

SNP_to_ensembl_one <- mapped_ld_snp_genes_csv %>%
  dplyr::select(SNP = SNP1, GENE_NAME = GENE_STABLE_1)

SNP_to_ensembl_combined <- rbind(SNP_to_ensembl_one, SNP_to_ensembl_second)
SNP_to_ensembl_combined_unique <- unique(SNP_to_ensembl_combined)

#all_linked_SNPs_DATA_genename <- all_linked_SNPs_DATA %>% dplyr::left_join(SNP_to_ensembl_combined_unique, by="SNP")
#all_UNLINKED_SNPs_DATA_genename <- all_UNLINKED_SNPs_DATA %>% dplyr::left_join(SNP_to_ensembl_combined_unique, by="SNP")
test <- all_linked_SNPs_DATA %>% dplyr::left_join(SNP_to_ensembl_combined_unique, by="SNP")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

gene_ids <- unique(SNP_to_ensembl_combined_unique$GENE_NAME) 
ens <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tss_tbl <- getBM(
  attributes = c("ensembl_gene_id",
                 "chromosome_name",
                 "strand",
                 "transcript_start",   # 5′ end of each transcript
                 "transcript_end"),    # 3′ end
  filters    = "ensembl_gene_id",
  values     = gene_ids,
  mart       = ens
) %>% 
  dplyr::mutate(
    tss_pos = ifelse(strand == 1, transcript_start, transcript_end)  # flip for minus strand
  ) %>% 
  dplyr::group_by(ensembl_gene_id, chromosome_name, strand) %>%            # multiple isoforms
  dplyr::summarise(tss_pos = min(tss_pos), .groups = "drop")   


snps  <- GRanges(seqnames = sig_snps_high_ES$CHR,
                 ranges   = IRanges(start = sig_snps_high_ES$POS, end = sig_snps_high_ES$POS),
                 SNP      = sig_snps_high_ES$SNP)

genes <- GRanges(seqnames = tss_tbl$chromosome_name,
                 ranges   = IRanges(start = tss_tbl$tss_pos, end = tss_tbl$tss_pos),
                 ensembl  = tss_tbl$ensembl_gene_id)

hits <- distanceToNearest(snps, genes)

nearest <- data.frame(
  SNP          = mcols(snps)$SNP[ queryHits(hits) ],
  ensembl_gene = mcols(genes)$ensembl[ subjectHits(hits) ],
  dist_bp      = mcols(hits)$distance
)

nearest_SNP_to_ensembl <- nearest %>% 
  dplyr::select(
    SNP, ensembl_gene
  )

all_linked_SNPs_DATA_genename <- all_linked_SNPs_DATA %>% dplyr::left_join(nearest_SNP_to_ensembl, by="SNP")
all_UNLINKED_SNPs_DATA_genename <- all_UNLINKED_SNPs_DATA %>% dplyr::left_join(nearest_SNP_to_ensembl, by="SNP")

all_SNPs_DATA_genename <- rbind(all_UNLINKED_SNPs_DATA_genename, all_linked_SNPs_DATA_genename)


#now let us look at SNP data
gene_vec <- unique(all_UNLINKED_SNPs_DATA_genename$ensembl_gene)
id_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys     = gene_vec,
                                keytype  = "ENSEMBL",
                                columns  = "ENTREZID")

entrez_genes <- unique(na.omit(id_map$ENTREZID))
length(entrez_genes)

universe_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys    = unique(all_SNPs_DATA_genename$ensembl_gene),
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID") %>%
  pull(ENTREZID) %>% na.omit() %>% unique()
length(universe_ids)

ego <- enrichGO(gene         = entrez_genes,
                #universe     = universe_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                minGSSize     = 3)   

n_over <- ego@result$geneID |>  
  str_split("/") |>   
  unlist() |>   
  unique() |> 
  length() 

unlinked <- dotplot(ego, showCategory = 11, font.size = 10) +
  ggtitle(paste0("Unlinked GO BP enrichment: ES > 0.010 (n = ", n_over, ")"))

ekegg <- enrichKEGG(
  gene          = entrez_genes,
  # universe   = universe_ids,
  organism      = "hsa",
  keyType       = "ncbi-geneid",  
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 2
)
n_over_kegg <- ekegg@result$geneID |>
  str_split("/") |>
  unlist() |>
  unique() |>
  length()

unlinked_kegg <- dotplot(ekegg, showCategory = 11, font.size = 10) +
  ggtitle(paste0("Unlinked KEGG enrichment: ES > 0.010 (n = ", n_over_kegg, ")"))

ereact <- enrichPathway(gene = entrez_genes, pvalueCutoff = 1, readable = TRUE)

# count unique genes in the enriched terms
n_over_react <- ereact@result$geneID |>
  strsplit("/") |>
  unlist() |>
  unique() |>
  length()

unlinked_react <- dotplot(ereact, showCategory = 20) +
  ggtitle(paste0("Reactome pathways (n = ", n_over_react, ")"))

#print(unlinked_react)
print(unlinked_kegg)

print(unlinked)

#
gene_vec <- unique(all_linked_SNPs_DATA_genename$ensembl_gene)
id_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys     = gene_vec,
                                keytype  = "ENSEMBL",
                                columns  = "ENTREZID")

entrez_genes <- unique(na.omit(id_map$ENTREZID))
length(entrez_genes)

universe_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys    = unique(all_SNPs_DATA_genename$ensembl_gene),
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID") %>%
  pull(ENTREZID) %>% na.omit() %>% unique()
length(universe_ids)

ego <- enrichGO(gene         = entrez_genes,
                #universe     = universe_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                minGSSize     = 3)     

n_over <- ego@result$geneID |>  
  str_split("/") |>   
  unlist() |>   
  unique() |> 
  length() 

linked <- dotplot(ego, showCategory = 11, font.size = 10) +
  ggtitle(paste0("Linked GO BP enrichment: ES > 0.010 (n = ", n_over, ")"))

ekegg <- enrichKEGG(
  gene          = entrez_genes,
  # universe   = universe_ids,   
  organism      = "hsa",        
  keyType       = "ncbi-geneid",  
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 2
)
n_over_kegg <- ekegg@result$geneID |>
  str_split("/") |>
  unlist() |>
  unique() |>
  length()

linked_kegg <- dotplot(ekegg, showCategory = 11, font.size = 10) +
  ggtitle(paste0("Linked KEGG enrichment: ES > 0.010 (n = ", n_over_kegg, ")"))

ereact <- enrichPathway(gene = entrez_genes, pvalueCutoff = 0.1, readable = TRUE)

# count unique genes in the enriched terms
n_over_react <- ereact@result$geneID |>
  strsplit("/") |>
  unlist() |>
  unique() |>
  length()

linked_react <- dotplot(ereact, showCategory = 20) +
  ggtitle(paste0("Reactome pathways (n = ", n_over_react, ")"))

#print(linked_react)

print(linked_kegg)

print(linked)

library(pathview) 
library(regioneR) 

kp_snps <- all_linked_SNPs_DATA_genename |>
  dplyr::mutate(CHR = paste0("chr", CHR))
chr_with_snps <- unique(kp_snps$CHR)
#chr_with_snps <- sortChromosomes(chr_with_snps)

kp <- plotKaryotype(genome = "hg38", chromosomes = chr_with_snps, plot.type = 2) 

maxLP <- ceiling(max(kp_snps$LP))
#kpAxis(kp, ymin = 0, ymax = maxLP, tick.pos = 0:maxLP, cex = 0.6)
kpAxis(kp,
       ymin = 0, ymax = maxLP) 

kpPoints(kp,
         chr  = kp_snps$CHR,
         x    = kp_snps$POS,
         y    = kp_snps$LP,
         ymin = 0, ymax = maxLP,
         col  = ifelse(kp_snps$P< 5e-8, "red", "steelblue"),
         cex  = 1.2, pch = 16)
kpAddMainTitle(kp, "Genome wide significant SNPs on HG38")

#UNLINKED
kp_snps <- all_UNLINKED_SNPs_DATA_genename |>
  dplyr::mutate(CHR = paste0("chr", CHR))
chr_with_snps <- unique(kp_snps$CHR)

kp <- plotKaryotype(genome = "hg38", chromosomes = chr_with_snps, plot.type = 4) 

maxLP <- ceiling(max(kp_snps$LP))
kpAxis(kp, ymin = 0, ymax = maxLP) 

kpPoints(kp,
         chr  = kp_snps$CHR,
         x    = kp_snps$POS,
         y    = kp_snps$LP,
         ymin = 0, ymax = maxLP,
         col  = ifelse(kp_snps$P< 5e-8, "red", "steelblue"),
         cex  = 1.2, pch = 16)
kpAddMainTitle(kp, "Genome wide significant UNLINKED SNPs on HG38")

#LINKED
kp_snps <- all_linked_SNPs_DATA_genename |>
  dplyr::mutate(CHR = paste0("chr", CHR))

kp <- plotKaryotype(genome = "hg38", chromosomes = chr_with_snps, plot.type = 4) 

maxLP <- ceiling(max(kp_snps$LP))
kpAxis(kp,
       ymin = 0, ymax = maxLP) 

kpPoints(kp,
         chr  = kp_snps$CHR,
         x    = kp_snps$POS,
         y    = kp_snps$LP,
         ymin = 0, ymax = maxLP,
         col  = ifelse(kp_snps$P< 5e-8, "red", "steelblue"),
         cex  = 1.2, pch = 16)
kpAddMainTitle(kp, "Genome wide significant LINKED SNPs on HG38")

library(STRINGdb)
ensembl_ids <- unique(all_UNLINKED_SNPs_DATA_genename$ensembl_gene)
gost_res <- gost(query     = ensembl_ids,
                 organism  = "hsapiens",
                 correction_method = "fdr")


head(gost_res$result[, c("term_name", "p_value", "source")])
gp_unlinked <- gostplot(gost_res, capped = TRUE)
print(gp_unlinked)




ensembl_ids <- unique(all_linked_SNPs_DATA_genename$ensembl_gene)
gost_res_linked <- gost(query     = ensembl_ids,
                 organism  = "hsapiens",
                 correction_method = "fdr")

# View top enriched terms
head(gost_res_linked$result[, c("term_name", "p_value", "source")])
gp_linked <- gostplot(gost_res_linked, capped = TRUE)
print(gp_linked)

N_TOP <- 3 
tbl <- gost_res_linked$result 

tbl_top <- tbl %>%
  filter(p_value < 0.05) %>%  
  group_by(source) %>%
  slice_max(order_by = -log10(p_value), n = N_TOP, with_ties = FALSE) %>%
  ungroup()

source_order <- c("GO:BP", "GO:MF", "GO:CC",
                  "KEGG", "REAC", "CORUM", "TF",
                  "MIRNA", "HP", "WP")
tbl_top <- tbl_top %>%
  mutate(source     = factor(source, levels = source_order),
         term_name  = str_trunc(term_name, 45)) 

p_bar <- tbl_top %>%
  mutate(term_name = fct_reorder(term_name, -log10(p_value))) %>%
  ggplot(aes(x = term_name,
             y = -log10(p_value),
             fill = source)) +
  geom_col(width = 0.7) +
  coord_flip() +
  geom_text(aes(label = round(-log10(p_value), 1)),
            hjust = -0.1, size = 3.2) +                      # value labels
  scale_fill_brewer(palette = "Set2", name = "") +
  labs(x = NULL,
       y = expression(-log[10]~"(adj. p)"),
       title = "Enrichment of Linked Genes (g:Profiler)") +
  theme_linedraw(base_size = 14) +
  theme(legend.position = "right",
        plot.title      = element_text(face = "bold"),
        axis.text.y     = element_text(size = 10))



gene_vec <- unique(all_SNPs_DATA_genename$ensembl_gene)
id_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys     = gene_vec,
                                keytype  = "ENSEMBL",
                                columns  = "ENTREZID")

entrez_genes <- unique(na.omit(id_map$ENTREZID))
length(entrez_genes)

universe_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys    = unique(all_SNPs_DATA_genename$ensembl_gene),
                                      keytype = "ENSEMBL",
                                      columns = "ENTREZID") %>%
  pull(ENTREZID) %>% na.omit() %>% unique()
length(universe_ids)

ego <- enrichGO(gene         = entrez_genes,
                #universe     = universe_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                minGSSize     = 3)     

n_over <- ego@result$geneID |>  
  str_split("/") |>   
  unlist() |>   
  unique() |> 
  length() 

all<- dotplot(ego, showCategory = 11, font.size = 10) +
  ggtitle(paste0("Linked/Unlinked GO BP enrichment: ES > 0.010 (n = ", n_over, ")"))
print(all)

