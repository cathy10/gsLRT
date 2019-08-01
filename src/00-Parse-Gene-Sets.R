## Load libraries
library(plyr)
library(tidyverse)
library(magrittr)
library(MASS)
## Read protein expression data
prot_data <- readr::read_csv("data/ADNI AD vs CN new.csv")
## Read gene sets
gs_list <- c("data/h.all.v6.2.symbols.gmt",
             "data/c2.cp.biocarta.v6.2.symbols.gmt",
             "data/c2.cp.kegg.v6.2.symbols.gmt",
             "data/c2.cp.reactome.v6.2.symbols.gmt",
             "data/c2.cp.v6.2.symbols.gmt",
             "data/c5.all.v6.2.symbols.gmt") %>%
  plyr::alply(1, qusage::read.gmt) %>%
  do.call(c, .)
## Read gene aliases
gene_alias <- readr::read_csv("data/gene_aliasing.csv")
## Create and rename protein expression matrix
G <- dplyr::select(prot_data, -RID, -sex, -age, -APOE, -Diagnosis) %>%
  as.matrix() %>%
  magrittr::set_colnames(gsub("_HUMAN", "", colnames(.)))
alias_colnames <- plyr::aaply(1:length(colnames(G)), 1, function(j){
  
  if(colnames(G)[j] %in% gene_alias[["gene_name"]]){
    gene_alias[["alias"]][gene_alias[["gene_name"]] == colnames(G)[j]]
  } else {
    colnames(G)[j]
  }
  
})
G %<>% magrittr::set_colnames(alias_colnames)
saveRDS(object = G, file = "results/objects/G.rds")
## Create disease status indicator
Z <- dplyr::select(prot_data, Diagnosis) %>%
  magrittr::extract2("Diagnosis") %>%
  magrittr::equals("AD") %>%
  ifelse(1, 0)
saveRDS(object = Z, file = "results/objects/Z.rds")
## Create matrix of additional covariates
X <- dplyr::select(prot_data, sex, age, APOE) %>%
  dplyr::mutate(sex = ifelse(sex == "F", 1, 0),
                age = as.numeric(gsub(" ", "", age)),
                APOE = ifelse(APOE == "E4", 1, 0)) %>%
  as.matrix()
saveRDS(object = X, file = "results/objects/X.rds")
## Trim gene set list
gene_names <- colnames(G)
gs_list %<>% names() %>%
  plyr::alply(1, function(n){
    gs <- gs_list[[n]]
    ind <- gs %in% gene_names
    if(sum(ind) < 5 || sum(ind) >= 100){return(NULL)}else{return(sort(gs[ind]))}
  }) %>%
  magrittr::set_names(names(gs_list))
gs_list[sapply(gs_list, is.null)] <- NULL
gs_list <- unique(gs_list) %>% 
  magrittr::set_names(names(gs_list)[!duplicated(gs_list)]) %>%
  magrittr::set_names(., gsub(".*[.]", "", names(.)))
saveRDS(object = gs_list, file = "results/objects/gene_sets.rds")
