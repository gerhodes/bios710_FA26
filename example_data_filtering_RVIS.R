library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)

my_working_dir <- ""

setwd(my_working_dir)

source("example_dat_filtering_helpers.R")

## read in current CCDS file, downloaded from: https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
ccds <- read.delim(
  "CCDS.current.txt",
  comment.char = "#",
  header = FALSE
)
## match column names to documentation
names(ccds) <- c("chromosome", "nc_accession", "gene", "gene_id", "ccds_id", "ccds_status",
                 "cds_strand", "cds_from", "cds_to", "cds_locations", "match_type")

## restrict to "Public" status on subset chromosomes
subset_chr <- c(19, 20, 21)
ccds_subset <- subset(ccds, (chromosome %in% as.character(subset_chr))&(ccds_status=="Public"))
head(ccds_subset$cds_locations)

## test parser
parse_cds(ccds_subset$cds_locations[1])

## convert all the transcripts into a GRanges object --> more memory efficient
gr_list <- lapply(
  seq_len(nrow(ccds_subset)),
  function(i) {
    
    coords <- parse_cds(ccds_subset$cds_locations[i])
    
    GRanges(
      seqnames = paste0("chr", ccds_subset$chromosome[i]),
      ranges = IRanges(
        start = coords$start,
        end   = coords$end
      ),
      gene = ccds_subset$gene[i],
      ccds_id = ccds_subset$ccds_id[i]
    )
  }
)

gr <- do.call(c, gr_list)

table(seqnames(gr))

## extend start and end positions by 2 bp to handle splice variants
start(gr) <- start(gr) - 2
end(gr) <- end(gr) + 2

## convert to a GRangesList by splitting on HGNC symbol
gene_gr <- split(gr, gr$gene)
gene_gr[["POTED"]] ## look at one example gene

## take the union of all CCDS transcripts per gene
gene_gr <- endoapply(gene_gr, reduce)
## make sure we still have intervals on all chromosomes
table(seqnames(unlist(gene_gr)))

## compute the coding sequence length for each gene
gene_lengths <- sapply(
  gene_gr,
  function(x) sum(width(x))
)

## read in coverage data, processed in separate script 
coverage <- readRDS("gnomad.v4.1.chr19_20_21_near_complete_callable.rds")

## convert to a GRanges object
coverage_gr <- GRanges(seqnames = coverage$chromosome,
                       ranges = IRanges(start=as.numeric(coverage$pos),
                                        end = as.numeric(coverage$pos)))
coverage_gr <- reduce(coverage_gr)

## compute coverage per gene
coverage_summ <- data.frame(gene = names(gene_gr),
                            coding_bp = gene_lengths,
                            covered_bp = numeric(length(gene_gr)))
## for each gene, find the overlap between the CDS and the near-complete callability regions
for(i in seq_along(gene_gr)) {
  overlaps <- intersect(
    gene_gr[[i]],
    coverage_gr
  )
  coverage_summ$covered_bp[i] <-
    sum(width(overlaps))
}
head(coverage_summ)
## compute the fraction covered
coverage_summ$fraction_covered <- coverage_summ$covered_bp / coverage_summ$coding_bp
summary(coverage_summ)

## define gene as "accessable" if at least 70% is covered
coverage_summ$accessable_flag <- coverage_summ$fraction_covered >= 0.7
head(coverage_summ)

## get vector with accessable genes
accessable_genes <- coverage_summ$gene[coverage_summ$accessable_flag]

## read in gnomAD chr19-21 variants
## note: this is a parsed version of the exome variant sites VCF file from the gnomAD browser
## which has already filtered out variants that do not pass QC filtering and to canonical transcript only
variants <- rbind(read.delim("gnomad_chr19_exome_sites_parsed.tsv", stringsAsFactors = F),
                  read.delim("gnomad_chr20_exome_sites_parsed.tsv", stringsAsFactors = F),
                  read.delim("gnomad_chr21_exome_sites_parsed.tsv", stringsAsFactors = F))

## filter to SNVs, remove large dataframe for efficiency
## SNV: both reference and alt allele will be 1 bp long
variants_snvs <- variants %>% filter(nchar(REF_x)==1, nchar(ALT_x)==1)
rm(variants)

## turn gnomAD variant data into a GRanges object
variants_gr <- GRanges(seqnames = variants_snvs$CHROM_x,
                       ranges = IRanges(start = variants_snvs$POS_x,
                                        end = variants_snvs$POS_x))

## subset CCDS to accessable genes
accessable_gene_gr <- gene_gr[intersect(names(gene_gr), accessable_genes)]

## convert back to one GRanges, retain gene names
gene_regions <- unlist(accessable_gene_gr)
gene_regions$gene <- rep(names(accessable_gene_gr), lengths(accessable_gene_gr))

## find variants in accessable gene CDS
overlap <- findOverlaps(variants_gr, gene_regions, ignore.strand=T)
queryHits(overlap)
subjectHits(overlap)

## generate variant-to-gene table
variant_gene_df <- data.frame(variant_idx = queryHits(overlap),
                              gene = gene_regions$gene[subjectHits(overlap)])

variants_in_accessable_genes <- cbind(variants_snvs[variant_gene_df$variant_idx, ],
                                      gene = variant_gene_df$gene)

## check 1: how many genes have at least one variant ?
length(unique(variants_in_accessable_genes$gene))

## check 2: are the VEP consequences we see what we expect for coding sequence ?
table(variants_in_accessable_genes$CONSEQUENCE)

## re-create the RVIS paper's functional vs non-functional classification
functional_consq <- c("missense_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant")
nonfunctional_consq <- c("synonymous_variant")

variants_in_accessable_genes <- variants_in_accessable_genes %>%
  mutate(
    classification = case_when(sapply(strsplit(CONSEQUENCE, "&"), function(x) any(x %in% functional_consq)) ~ "functional",
                               sapply(strsplit(CONSEQUENCE, "&"), function(x) any(x %in% nonfunctional_consq)) ~ "non_functional",
                               TRUE ~ NA_character_)
  )

## check classification worked properly
variants_in_accessable_genes %>% filter(classification=="functional") %>% count(CONSEQUENCE)
variants_in_accessable_genes %>% filter(classification=="non_functional") %>% count(CONSEQUENCE)
variants_in_accessable_genes %>% filter(is.na(classification)) %>% count(CONSEQUENCE)

## keep only variants with a classification 
variants_in_accessable_genes <- variants_in_accessable_genes %>% filter(!is.na(classification))

## generate gene counts for regression
## Y: number of common functional variants, define common to be MAF>0.1%
## X: total number of coding variants, functional or non-functional regardless of MAF

gene_counts <- variants_in_accessable_genes %>% group_by(gene) %>%
  summarise(
    num_common_functional = sum((AF_all>0.001)&(classification=="functional")),
    total_coding_variants = n()
  )

gene_counts <- gene_counts %>% filter(total_coding_variants < 10000)

ggplot(data=gene_counts, mapping=aes(x=total_coding_variants, y=num_common_functional)) +
  geom_point() +
  xlab("Total Num. Protein-Coding Variants") +
  ylab("Num. Common (MAF>0.1%) Functional Variants") +
  theme_bw() +
  ggtitle("Total Protein-Coding Variants vs Common Functional Variants: chr19, chr20, chr21")

saveRDS(gene_counts, "chr19_20_21_gnomADv4_gene_counts_for_RVIS.rds")
write_tsv(gene_counts, "gnomAD_gene_counts_subset_RVIS.tsv")
