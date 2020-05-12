library(Biostrings)
library(dplyr)
d <- Biostrings::readDNAStringSet("ncbi_sequences.fasta")

ssd <- strsplit(names(d), split="|", fixed=T)
meta_data <- data.frame(sequence_id=sapply(ssd, "[", 1),
                        country=sapply(ssd, "[", 4),
                        length=nchar(d),
                        index=1:length(d),
                        stringsAsFactors = F) %>%
            dplyr::filter(!grepl(" ", country), length > 20000)

# downsample the USA genomes
USA_ind <- which(meta_data$country == "USA")
USA_ind_sample <- USA_ind[sample.int(length(USA_ind), 100)]

meta_data <- rbind(meta_data[USA_ind_sample,],
                   dplyr::filter(meta_data, country != "USA"))


d_filtered <- d[meta_data$ind]
ref <- df[1] %>% as.character

all_patterns <- rep("", length(d_filtered))
for (i in 1:length(d_filtered)) {
  cat("alignment", i, "of", length(d_filtered), "\n")
  pa <- pairwiseAlignment(pattern = as.character(d_filtered[i]), subject=ref)
  pattern_a <- strsplit(alignedPattern(pa) %>% as.character, split="")[[1]]
  subject_a <- strsplit(alignedSubject(pa) %>% as.character, split="")[[1]]

  aligned_pattern <-  pattern_a[subject_a != "-"]
  all_patterns[i] <- paste(aligned_pattern, collapse="")
}

out_df <- DNAStringSet(all_patterns)
names(out_df)  <-paste("sequence", 1:length(d_filtered), "_", meta_data$country, sep="")
writeXStringSet(out_df, "ncbi_sequences_aligned.fasta")




