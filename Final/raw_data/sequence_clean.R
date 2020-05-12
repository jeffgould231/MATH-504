library(Biostrings); library(dplyr); library(data.table)
dr <- readDNAStringSet("reference.fasta")
d <- readDNAMultipleAlignment("ncbi_sequences_aligned.fasta")
m <- as.matrix(d)

countries <- sapply(strsplit(rownames(m), split="_"), "[", 2)

# remove N and -
ref <- strsplit(as.character(dr), split="")[[1]]
m_clean <- apply(m, 1, function(r) {
  bad_ind <- r=="N" | r=="-"
  r[bad_ind] <- ref[bad_ind]
  return (r)
}) %>% t


fwrite(m_clean, "matrix_sequences.csv", col.names=F)
fwrite(matrix(countries, ncol=1), "matrix_countries.csv", col.names=F, quote=T)
fwrite(matrix(ref, nrow=1), "matrix_reference.csv", col.names=F)

final_d <- DNAStringSet(apply(m_clean, 1, paste, collapse=""))
names(final_d) <- paste("sequence", 1:length(final_d), "_", countries, sep="")
writeXStringSet(final_d, "ncbi_sequences_cleaned.fasta")
