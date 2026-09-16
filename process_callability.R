library(tidyverse)

my_working_dir <- ""

setwd(my_working_dir)

## read in all sites allele numbers, downloaded from: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/tsv/exomes/gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz
coverage <- read_tsv("gnomad.exomes.v4.1.allele_number_all_sites.tsv.bgz")

## split locus column on ":" to get chromosome and genomic location columns
coverage <- coverage %>% separate_wider_delim(cols=locus,
                                              delim=":",
                                              names=c("chromosome", "pos"))

## filter to subset chromosomes, remove large all sites dataframe from environment
# subset_chr <- c(19, 20, 21)
# subset_chr <- seq(1:18)
subset_chr <- 22
coverage_subset <- coverage %>% filter(chromosome %in% paste0("chr", subset_chr))
rm(coverage)

## determine the maximum AN
max_AN <- max(coverage_subset$AN, na.rm = TRUE)

## define a site as "callable" if the AN is at least 90% of the maximum AN
coverage_subset_callable <- coverage_subset %>% filter(AN >= 0.9*max_AN)

## save the callable sites df
saveRDS(coverage_subset_callable, "gnomad.v4.1.chr22_near_complete_callable.rds")
