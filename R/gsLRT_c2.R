require(magrittr)
require(plyr)
require(tidyverse)
source('gsLRT.R')
source('lr_stats.R')
source('plot_results.R')
source('chisq_test.R')
source('enrich_score.R')

## Read R objects
G <- readRDS("G_jordan.rds")
Z <- readRDS("Z_jordan.rds")
X <- readRDS("X_jordan.rds") #good
gene_sets <- readRDS("gene_sets_jordan.rds") 

#running gsLRT function
ES_df <- gsLRT(Z, X, G, interaction_term = "sex", gene_sets = gene_sets, n_perm = 200)

########
G1 <- readRDS("G.rds")
Z1 <- readRDS("Z.rds")
X1 <- readRDS("X_sub.rds")
gene_sets1 <- readRDS("gene_sets.rds") 

# setting up my X
X2 <- X1[c(3,2,4,5)]
colnames(X2) <- c('Condition', 'HN', 'sex', 'age') #final X value

X2_1 <- X2[2:4] #taking out the Condition column
Matrix_X <- matrix(unlist(X2_1), ncol = 3, byrow = F) #converting the list to a matrix
colnames(Matrix_X) <- c('HN', 'sex', 'age')  #renaming them

# setting up my G
Matrix_G <- matrix(unlist(G1), ncol = 2014, byrow = F) #converting the list to a matrix
colnames(Matrix_G) <- gene_names  #renaming them

# setting up my Z
Matrix_Z <- as.numeric(unlist(Z1))

# setting up my gene_sets
Matrix_gene_sets <- lapply(gene_sets_c2, as.character) 
  
  
#running gsLRT function
ES_df <- gsLRT(Matrix_Z, Matrix_X, Matrix_G, interaction_term = "sex", gene_sets = Matrix_gene_sets, n_perm = 200)





