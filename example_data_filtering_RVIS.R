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
# subset_chr <- c(19, 20, 21)
subset_chr <- seq(1:22)
# subset_chr <- 22
ccds_subset <- subset(ccds, (chromosome %in% as.character(subset_chr))&(ccds_status=="Public"))
table(ccds_subset$chromosome)
head(ccds_subset$cds_locations)

## test parser
parse_cds(ccds_subset$cds_locations[1])

## convert all the transcripts into a GRanges object --> more memory efficient
gr_list <- lapply(
  seq_len(nrow(ccds_subset)),
  function(i) {
    
    coords <- parse_cds(ccds_subset$cds_locations[i])
    if(anyNA(coords)){
      print("Parsing produced NAs:")
      print(i)
      print(ccds_subset$cds_locations[i])
      print(coords)
      
    }
    if(!anyNA(coords)){
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
    
  }
)

gr <- do.call(c, gr_list)

table(seqnames(gr))
table(ccds_subset$chromosome)

## extend start and end positions by 2 bp to handle splice variants
start(gr) <- start(gr) - 2
end(gr) <- end(gr) + 2

## convert to a GRangesList by splitting on HGNC symbol
gene_gr <- split(gr, gr$gene)
gene_gr[["BRCA1"]] ## look at one example gene

## take the union of all CCDS transcripts per gene
gene_gr <- endoapply(
  gene_gr,
  GenomicRanges::reduce
)

## make sure we still have intervals on all chromosomes
table(seqnames(unlist(gene_gr)))

## compute the coding sequence length for each gene
gene_lengths <- sapply(
  gene_gr,
  function(x) sum(width(x))
)

## read in coverage data, processed in separate script 
coverage <- rbind(readRDS("gnomad.v4.1.chr1_to_18_near_complete_callable.rds"),
                  readRDS("gnomad.v4.1.chr19_20_21_near_complete_callable.rds"),
                  readRDS("gnomad.v4.1.chr22_near_complete_callable.rds"))
# coverage <- readRDS("gnomad.v4.1.chr22_near_complete_callable.rds")

## convert to a GRanges object
coverage_gr <- GRanges(seqnames = coverage$chromosome,
                       ranges = IRanges(start=as.numeric(coverage$pos),
                                        end = as.numeric(coverage$pos)))
coverage_gr <- GenomicRanges::reduce(coverage_gr)

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
sum(coverage_summ$accessable_flag) ## num. accessable genes

## get vector with accessable genes
accessable_genes <- coverage_summ$gene[coverage_summ$accessable_flag]

## read in gnomAD variant data (SNVs classified across transcripts)
variants_unique <- readRDS("gnomad_v4.1_unique_classified_SNVs.rds")

## turn gnomAD variant data into a GRanges object
variants_gr <- GRanges(seqnames = variants_unique$CHROM_x,
                       ranges = IRanges(start = variants_unique$POS_x,
                                        end = variants_unique$POS_x))

## subset CCDS to accessable genes
accessable_gene_gr <- gene_gr[accessable_genes]

## convert back to one GRanges, retain gene names
gene_regions <- unlist(accessable_gene_gr)
gene_regions$gene <- rep(names(accessable_gene_gr), lengths(accessable_gene_gr))

## find variants in accessable gene CDS
overlap <- findOverlaps(variants_gr, gene_regions, ignore.strand=T)
queryHits(overlap)
subjectHits(overlap)

## generate variant-to-gene table
variant_gene_df <- data.frame(variant_id = variants_unique$variant_id[queryHits(overlap)],
                              gene = gene_regions$gene[subjectHits(overlap)])


variants_final <- variant_gene_df %>%
  left_join(
    variants_unique,
    by = "variant_id"
  )

## check 1: how many genes have at least one variant ?
length(unique(variants_final$gene))

saveRDS(variants_final, "gnomAD_variants_CDS_acessable_classified.rds")

## generate gene counts for regression
## Y: number of common functional variants, define common to be MAF>0.1%
## X: total number of coding variants, functional or non-functional regardless of MAF

# variants_final <- readRDS("gnomAD_variants_CDS_acessable_classified.rds")

gene_counts <- variants_final %>% group_by(gene) %>%
  summarise(
    num_common_functional = sum((AF_all>0.001)&(classification=="functional")),
    total_coding_variants = n()
  )

gene_counts_filt <- gene_counts %>% filter(num_common_functional < 140)

ggplot(data=gene_counts_filt, mapping=aes(x=total_coding_variants, y=num_common_functional)) +
  geom_point() +
  xlab("Total Num. Protein-Coding Variants") +
  ylab("Num. Common (MAF>0.1%) Functional Variants") +
  theme_bw() +
  ggtitle("Total Protein-Coding Variants vs Common Functional Variants: gnomAD v4")

write_tsv(gene_counts, "gnomAD_gene_counts_subset_RVIS_upd.tsv")
