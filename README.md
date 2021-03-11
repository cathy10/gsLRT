# gsLRT
Gene set enrichment from likelihood ratio statistics

# Running gsLRT on small scale genomics data

To run gsLRT on your own genomics data, open the R console or an RStudio window and follow the steps below. To make the code included here compatible with your data with only minor modifications you will need the following:

- An `.rds` object `G` containing a matrix of values from some genomics assay (e.g. gene expression, protein expression). The columns of this matrix should be labeled with gene symbols (symbols must match those in the gene sets you use later) and the rows of this matrix should correspond to your experimental samples. For instance, for a an assay measuring expression of 25 genes collected from 10 mice, this object would be a 10 x 25 matrix.
- An `.rds` object `Z` containing a binary label for each sample. Following the example presented above, this should be a length 10 vector of zeros and ones, indicating absence/presence of some phenotype of interest.
- An `.rds` object `X` containing a matrix of `k` covariates for each sample. For instance, a 10 x 3 matrix containing age, sex (i.e. 0 for male, 1 for female), and Mut/Wt status for a gene (also a 0, 1 variable).
- An `.rds` object `gene_sets` containing a list of gene sets. Each element of the list is named by a gene set name, e.g. "GO_STRUCTURAL_MOLECULE_ACTIVITY", and each element is itself a vector of strings, where each string is a unique gene symbol.

With these objects in hand, here is an example `gsLRT` run. In this example, we are specifically testing for whether an interaction between the phenotype `Z` and sex is significantly associated with any gene sets.

First, read in the data:

```
require(magrittr)
require(plyr)
require(tidyverse)

## Read R objects
G <- readRDS("results/objects/G.rds")
Z <- readRDS("results/objects/Z.rds")
X <- readRDS("results/objects/X.rds")
gene_sets <- readRDS("results/objects/gene_sets.rds")
```

Next, run the permutation-based version of `gsLRT` with

```
ES_df <- gsLRT::gsLRT(Z, X, G, interaction_term = "sex", gene_sets = gene_sets, n_perm = 2000)
```

This gives you a `data.frame` of gene sets and their permutation p-values. While the permutation-based test is recommended, depending on the size of your data, it may be prohibitively slow. In that case, a less conservative chi-squared test may be used:

```
Lambdas <- plyr::aaply(1:ncol(G), 1, function(j){
  
  gsLRT::lr_stats(y = G[, j], X = X, z = Z, interaction_term = "sex",
                  stat = "bartlett")
  
}, .progress = "text") %>%
  t() %>%
  magrittr::set_colnames(colnames(G))
ES <- gsLRT::enrich_score(Lambdas, gene_sets, trim = T, 
                           min_size = 5, stat = "chisq")
ES_df <- gsLRT::chisq_test(ES, gene_sets)

## Plot a figure
gsLRT::plot_results(ES_df, Lambdas, gene_sets)
```


