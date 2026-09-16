library(tidyverse)
## read in gnomAD variants
## note: this is a parsed version of the exome variant sites VCF file from the gnomAD browser
## which has already filtered out variants that do not pass QC filtering and to canonical transcript only
variants <- rbind(read.delim("gnomad_data/chr1.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr2.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr3.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr4.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr5.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr6.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr7.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr8.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr9.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr10.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr11.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr12.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr13.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr14.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr15.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr16.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr17.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr18.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr19.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr20.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr21.final.tsv", stringsAsFactors = F),
                  read.delim("gnomad_data/chr22.final.tsv", stringsAsFactors = F))
# variants <- read.delim("gnomad_chr22_exome_sites_parsed.tsv", stringsAsFactors = F)

## filter to SNVs, remove large dataframe for efficiency
## SNV: both reference and alt allele will be 1 bp long
variants_snvs <- variants %>% filter(nchar(REF_x)==1, nchar(ALT_x)==1)
rm(variants)

## create unique variant IDs
variants_snvs <- variants_snvs %>% mutate(variant_id = paste(CHROM_x, POS_x, REF_x, ALT_x, sep=":"))
head(variants_snvs)

## re-create the RVIS paper's functional vs non-functional classification
functional_consq <- c("missense_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant")
nonfunctional_consq <- c("synonymous_variant")

## classify each individual transcript record
## note: gnomAD data may include multiple transcripts per variant, we will classify all transcripts as functional
## or non-functional and determine per variant-id whether the variant is functional anywhere
## using functional > non-functional > ignore
variants_snvs <- variants_snvs %>%
  mutate(
    transcript_class = case_when(sapply( strsplit(CONSEQUENCE, "&"),
                                         function(x)
                                           any(x %in% functional_consq)) ~ "functional",
                                 sapply(strsplit(CONSEQUENCE, "&"),
                                        function(x)
                                          any(x %in% nonfunctional_consq)) ~ "non_functional",
                                 TRUE ~ NA_character_))

## check classification worked properly
variants_snvs %>% filter(transcript_class=="functional") %>% count(CONSEQUENCE)
variants_snvs %>% filter(transcript_class=="non_functional") %>% count(CONSEQUENCE)
variants_snvs %>% filter(is.na(transcript_class)) %>% count(CONSEQUENCE)

## collapse to one classification per variant
## first check: if functional on any transcript, classify the variant as functional
## second check: if it isn't functional on any transcript, then if it is non-functional on any transcript classify as non-functional
## third check: if it is not functional or non-functional on any transcript, classify as NA
variants_snvs <- variants_snvs %>%
  mutate(
    class_rank = case_when(
      transcript_class == "functional" ~ 2L,
      transcript_class == "non_functional" ~ 1L,
      TRUE ~ 0L
    )
  )

variant_snvs_classes <- variants_snvs %>%
  group_by(variant_id) %>%
  summarise(
    class_rank = max(class_rank),
    .groups = "drop"
  ) %>%
  mutate(
    classification = case_when(
      class_rank == 2L ~ "functional",
      class_rank == 1L ~ "non_functional",
      TRUE ~ NA_character_
    )
  )
`
## reduce to one row per variant
variants_snvs_unique <- variants_snvs %>%
  distinct(variant_id, CHROM_x, POS_x, REF_x, ALT_x,
           AF_all, AC_all, AN_all)

## attach the classification
variants_snvs_unique <- variants_snvs_unique %>%
  left_join(
    variant_snvs_classes,
    by = "variant_id"
  )

## drop variants not classified as functional or non-functional
variants_snvs_unique <- variants_snvs_unique %>%
  filter(!is.na(classification))

## save unique variants file
saveRDS(variants_snvs_unique, "gnomad_v4.1_unique_classified_SNVS.rds")
