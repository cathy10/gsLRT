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
saveRDS(Matrix_G, "Matrix_G.rds") #needed to convert the mice genes to human 
#retrieve human protein names matrix G
Matrix_G_final = readRDS("Matrix_G_final.rds")

# setting up my Z
Matrix_Z <- as.numeric(unlist(Z1))

# setting up my gene_sets
gene_sets_c2 <- readRDS("gene_sets_c2.rds")

#change the gene_sets_c2.rds to a list of lists
Matrix_gene_sets <- vector(mode = "list", length = 40)
  ##get the list of pathway names
pathway_name = rep(NA, 40)
pathway_name = as.character(gene_sets_c2[[1]])
  ##assembling the pathway names and corresponding proteins 
for (i in 1:length(pathway_name)){
  Matrix_gene_sets[[i]] <- as.character(gene_sets_c2[i, 2:length(gene_sets_c2[[i]])])
}
names(Matrix_gene_sets) <- pathway_name #rename the lists with pathways
saveRDS(Matrix_gene_sets, "c2_gene_sets_beforeconvert.rds")
#need to take out the empty strings

#this one contains empty pathways
Matrix_gene_sets_converted <- readRDS("Matrix_c2gene_sets_converted.rds")
#this one contains empty pathways
Matrix_gene_sets_morethan2 <- readRDS("Matrix_gene_sets_morethan2.rds")


#manipulate Matrix_X to be the same as X_jordan.rds
##reorder the matrix columns
Matrix_X1 = Matrix_X[, c(2, 3, 1)]
colnames(Matrix_X1) = c("sex", "age", "APOE")

#save all the rds files ready to run the gsLRT
saveRDS(Matrix_X1, "Matrix_X1.rds")
saveRDS(Matrix_Z, "Matrix_Z.rds")
saveRDS(Matrix_G_final, "Matrix_G_final.rds")


#running gsLRT function
ES_df1 <- gsLRT(Matrix_Z, Matrix_X1, Matrix_G_final, interaction_term = "age", 
               gene_sets = Matrix_gene_sets_morethan2, n_perm = 200)




