require(magrittr)
require(plyr)
require(tidyverse)
source('gsLRT.R')
source('lr_stats.R')
source('plot_results.R')
source('chisq_test.R')
source('enrich_score.R')


G <- readRDS("Matrix_G_final.rds")
Z <- readRDS("Matrix_Z.rds")
X <- readRDS("Matrix_X1.rds")
gene_sets <- readRDS("Matrix_gene_sets_morethan2.rds") 

ES_df <- gsLRT(Z, X, G, interaction_term = "sex", gene_sets = gene_sets, n_perm = 200)
