#
require(magrittr)
require(plyr)
require(tidyverse)
# Run gsLRT with no covariates and no interaction terms
## Read R objects
G <- readRDS("results/objects/G.rds")
Z <- readRDS("results/objects/Z.rds")
X <- readRDS("results/objects/X.rds")
gene_sets <- readRDS("results/objects/gene_sets.rds")
## Run gsLRT
Lambdas <- plyr::aaply(1:ncol(G), 1, function(j){
  
  gsLRT::lr_stats(y = G[, j], X = X, z = Z, interaction_term = NULL,
                  stat = "bartlett")
  
}, .progress = "text") %>%
  t() %>%
  magrittr::set_colnames(colnames(G))
##
ES_df <- gsLRT::gsLRT(Z, X, G, gene_sets = gene_sets, n_perm = 2000)
## Make figures
gsLRT::plot_results(ES_df, Lambdas, gene_sets)
ggplot2::ggsave("results/figures/perm_hyp_cov.jpeg", height = 10, width = 5)
## Write tables
ES_df %>%
  plyr::ddply(.(GS), function(d){
    gene_sets[[d[["GS"]]]] %>%
      sort() %>%
      paste(collapse = ", ") %>%
      data.frame(d, genes = .)
  }) %>%
  plyr::arrange(p) %>%
  readr::write_tsv("results/tables/perm_hyp_cov.txt")
