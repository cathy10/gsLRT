#convert mice and human gene names

#installing the required pacakge:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
library(tidyverse)
require("biomaRt")

#testing:
#musGenes <- c("Hmmr", "Tlx3", "Cpeb4")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = x , mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
#testing:
#genes <- convertMouseGeneList(musGenes)

gene_names <- readRDS("gene_names.rds") 
human_gene_names <- convertMouseGeneList(gene_names)



# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                   values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  no_mouse_genes <- length(x)
  no_human_genes <- length(humanx)
  
  if(no_human_genes != no_mouse_genes){
    print("Some genes could not be translated!")
    genes_not_trans <- setdiff(x,genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  # Print all gene names that could not be translated and the number of genes that were not translated
  
  return(humanx)
}


genes <- convertMouseGeneList(gene_names)
###
library(stringi)

f = read.delim("uniprot.tab", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
converter <- f %>% 
  dplyr::select(1, 3)
colnames(converter) = c("mice_gene", "human_gene")
clean <- gsub('_HUMAN', ' ', converter['human_gene'][[1]])
converter['human_gene_clean'] = clean
#cleaned data for converting gene names from mice to human
converter <- converter %>% 
  dplyr::select(1, 3)

gene_names <- readRDS('gene_names.rds')

#regardless of lower or upper case
# human_genes = c()
# micegene_nf = c()
# counter_yes = 0
# counter_no = 0
# for (micegene in gene_names){
#   if  (tolower(micegene) %in% tolower(converter[[1]])) {
#     found_gene = stri_detect_fixed(tolower(converter[[1]]), tolower(micegene))
#     human_genes = c(human_genes, converter[[2]][found_gene])
#     counter_yes = counter_yes + 1
#     #print(paste0("mice gene found: ", human_genes))
#   } else {
#     micegene_nf = c(micegene_nf, micegene)
#     counter_no = counter_no + 1
#   }
# }
print(paste0("Number converted to human gene: ", length(human_genes)))
print(paste0("Number of mice gene not found: ", length(micegene_nf)))

library(tidyverse)

human_genes = c()
micegene_nf = c()
counter_yes = 0
counter_no = 0
OneToMore_matching = tibble()
OneToMore_matching = OneToMore_matching %>%
  add_column(mice = '',
             human = '')

OneToOne_matched = tibble()
OneToOne_matched = OneToOne_matched %>%
  add_column(mice = '',
             human = '')
for (micegene in gene_names){
  if  (micegene %in% converter[[1]]) {
    found_gene = stri_detect_fixed(converter[[1]], micegene)
    if (length(converter[[2]][found_gene]) > 1){
      OneToMore_matching = OneToMore_matching %>% add_row(mice = micegene,
                                     human = converter[[2]][found_gene]) 
      #print(paste0(micegene, ": ", converter[[2]][found_gene]))
    } else {
      OneToOne_matched = OneToOne_matched %>% add_row(mice = micegene,
                                     human = converter[[2]][found_gene]) 
      }
    human_genes = c(human_genes, converter[[2]][found_gene])
    counter_yes = counter_yes + 1
    #print(paste0("mice gene found: ", human_genes))
  } else {
    micegene_nf = c(micegene_nf, micegene)
    counter_no = counter_no + 1
  }
}
library("readxl")
#generate the file containing the mice genes that can be translated to more than 1 human proteins.
write_csv(OneToMore_matching, "OneToMore_matching.csv")
OneToMore_cleaned = read_excel('OneToMore_cleaned.xlsx')
matched_genes = rbind(OneToMore_cleaned, OneToOne_matched) #final file for mice genes that are matched

#select the rolls in matrix_G that has human protein mapping
Matrix_G <- readRDS("Matrix_G.rds")

matched_gene_names <- matched_genes[[1]]
Matrix_G_cleaning <- Matrix_G[, matched_gene_names]
colnames(Matrix_G_cleaning) <- matched_genes[[2]]


#select columns that have gene match



