#
require(magrittr)
require(plyr)
require(tidyverse)
# Run gsHyp with no covariates and no interaction terms
## Read R objects
G <- readRDS("results/objects/G.rds")
Z <- readRDS("results/objects/Z.rds")
X <- readRDS("results/objects/X.rds")
gene_sets <- readRDS("results/objects/gene_sets.rds")
## Run gsHyp
Lambdas <- plyr::aaply(1:ncol(G), 1, function(j){
  
  gsHyp::lr_stats(y = G[, j], X = X, z = Z, interaction_term = "age",
                  stat = "bartlett")
  
}, .progress = "text") %>%
  t() %>%
  magrittr::set_colnames(colnames(G))
##
ES_df <- gsHyp::gsHyp(Z, X, G, interaction_term = "age", gene_sets = gene_sets, n_perm = 2000)
## Make figures
gsHyp::plot_results(ES_df, Lambdas, gene_sets)
ggplot2::ggsave("results/figures/perm_hyp_age.jpeg", height = 10, width = 5)
## Write tables
ES_df %>%
  plyr::ddply(.(GS), function(d){
    gene_sets[[d[["GS"]]]] %>%
      sort() %>%
      paste(collapse = ", ") %>%
      data.frame(d, genes = .)
  }) %>%
  plyr::arrange(p) %>%
  readr::write_tsv("results/tables/perm_hyp_age.txt")
