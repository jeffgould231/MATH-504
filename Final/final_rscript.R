###########################
sequences_load <- read_csv("matrix_sequences.csv", col_names = FALSE, col_types = cols(.default = "c"))
countries <- read_csv("matrix_countries.csv", col_names = FALSE)
reference <- read_csv("matrix_reference.csv", col_names = FALSE, col_types = cols(.default = "c"))

cl <- parallel::makeCluster(10)
parallel::clusterExport(cl = cl, varlist = "reference")
sequences <- t(parallel::parApply(cl = cl, sequences_load, 1, function(x){as.numeric(x == reference[1,])}))
parallel::stopCluster(cl)

mu <- colMeans(sequences)
sequences_center <- t(apply(sequences, 1, function(x){x-mu}))

n <- ncol(sequences_center)
N <- nrow(sequences_center)

V <- matrix(runif(2 * n), nrow = n)
V[,1] = V[,1] / norm(V[,1])
V[,2] = V[,2] / norm(V[,2])

for (i in 1:1000) {
  V <- t(sequences_center) %*% (sequences_center %*% V)
  V <- qr.Q(qr(V))
}

PCs <- sequences_center %*% V

plot_data <- data.frame(Country = countries, PC1 = PCs[,1], PC2 = PCs[,2]) %>%
  rename(Country = X1) 

variance <- sum(apply(sequences_center, 2, function(x){sum(x^2)}))
lambda_1 <- norm(sequences_center %*% V[,1]) ^2 / (norm(V[,1])^2)
lambda_2 <- norm(sequences_center %*% V[,2]) ^2 / (norm(V[,2])^2)
var_captured <- (lambda_1 + lambda_2) / variance

library(RColorBrewer)

ii <- length(unique(countries$X1)) %>% as.numeric()
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,ii), col=sample(col_vector, ii))

ggplot(data = plot_data, aes(x = PC1, y = PC2, color = Country)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = col_vector) +
  labs(title = "Principle Component Analysis of SARS-Cov2 Virus",
       x = glue::glue("PC1 ({round(lambda_1 / variance *100,2)}% of the variance)"),
       y = glue::glue("PC2 ({round(lambda_2 / variance *100,2)}% of the variance)"))






###############################







