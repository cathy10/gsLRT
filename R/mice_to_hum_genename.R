#convert mice and human gene names

#installing the required pacakge:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
library(tidyverse)
library(stringi)
library("readxl")

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

print(paste0("Number converted to human gene: ", length(human_genes)))
print(paste0("Number of mice gene not found: ", length(micegene_nf)))


#generate the file containing the mice genes that can be translated to more than 1 human proteins.
write_csv(OneToMore_matching, "OneToMore_matching.csv")
OneToMore_cleaned = read_excel('OneToMore_cleaned.xlsx')
matched_genes = rbind(OneToMore_cleaned, OneToOne_matched) #final file for mice genes that are matched

#select the columns in matrix_G that has human protein mapping
Matrix_G <- readRDS("Matrix_G.rds")

matched_gene_names <- matched_genes[[1]]
matched_gene_names.df = as.data.frame(matched_gene_names)
write_csv(matched_gene_names.df, "matched_gene_names.csv")

Matrix_G_cleaning <- Matrix_G[, matched_gene_names]
colnames(Matrix_G_cleaning) <- matched_genes[[2]]

saveRDS(Matrix_G_cleaning, "Matrix_G_final.rds") #ready to be used for gsLRT


#append converter with OneToMore_cleaned
##get unique mice genes in converted
colnames(OneToMore_cleaned) <- c("mice_gene", "human_gene_clean")
full_converter = rbind(OneToMore_cleaned, converter)
full_converter = full_converter[!duplicated(full_converter$mice_gene), ]
saveRDS(full_converter, "final_converted.rds")

################
#converting mice gene to human protein in the c2's pathways
unconverted = 0
converted = 0
empty_pathways = rep()
pathway_lengths = rep()
Matrix_gene_sets_c2 = readRDS('c2_gene_sets_beforeconvert.rds')
Matrix_gene_sets_converted = vector(mode = "list", length = 40)
for (i in 1:length(Matrix_gene_sets_c2)){
  for (micegene in Matrix_gene_sets_c2[[i]]){
    if (micegene != ""){
      if (micegene %in% toupper(full_converter$mice_gene)){
        human = full_converter$human_gene_clean[toupper(full_converter$mice_gene) == micegene]
        Matrix_gene_sets_converted[[i]] = c(Matrix_gene_sets_converted[[i]], human)
        converted = converted + 1
      } else {
        unconverted = unconverted + 1
      }
    }
  }
  #get a list of pathway index if the pathway contains two or less genes
  #if (is.null(Matrix_gene_sets_converted[[i]])) {empty_pathways = c(empty_pathways, i)}
  if (length(Matrix_gene_sets_converted[[i]]) < 3) {
    empty_pathways = c(empty_pathways, i)
  }
  pathway_lengths = c(pathway_lengths, length(Matrix_gene_sets_converted[[i]]))
}

#a dataframe containing the pathway names and number of proteins each
# pathways_df = data.frame(pathyway = pathway_name, protein_num = pathway_lengths)
# write_csv(pathways_df, 'pathway_df.csv')

names(Matrix_gene_sets_converted) <- names(Matrix_gene_sets_c2)

#Taking out the empty ones and the ones less than 2 genes
Matrix_gene_sets_non_null = Matrix_gene_sets_converted[-empty_pathways] #delete the ones with empty pathways
pathway_to_delete = rep()
pathway_lengths1 = rep()
pathway_lengths2 = rep()
for (i in 1:length(Matrix_gene_sets_non_null)){ 
  if (length(Matrix_gene_sets_non_null[[i]]) < 2) {
      pathway_to_delete = c(pathway_to_delete, i)
  }
    pathway_lengths1 = c(pathway_lengths1, length(Matrix_gene_sets_non_null[[i]]))
}
Matrix_gene_sets_morethan2 = Matrix_gene_sets_non_null#[-pathway_to_delete]

#get the unconverted pathway matrix to be of the same length as Matrix_gene_sets_morethan2
Matrix_gene_sets_c2_nonnull = Matrix_gene_sets_c2[-empty_pathways]
Matrix_gene_sets_unconverted_morethan2 = Matrix_gene_sets_c2_nonnull#[-pathway_to_delete]
#write to csv Matrix_gene_sets_c2_nonnull.df
write.table(Matrix_gene_sets_c2_nonnull, file = "Matrix_gene_sets_c2_nonnull.csv")
Matrix_gene_sets_c2_nonnull.df <- do.call("rbind", sapply(Matrix_gene_sets_c2_nonnull, as.data.frame)) 
write.csv(Matrix_gene_sets_c2_nonnull.df, file = "Matrix_gene_sets_c2_nonnull.csv")
pathway_name_c2_nonnull = names(Matrix_gene_sets_c2_nonnull)
write.csv(pathway_name_c2_nonnull, "pathway_name_c2_nonnull.csv")
#calculate the percentage of genes successfully converted to human in each pathway:
#num of proteins converted found in the pathway/num of protein in the pathway total before conversion
percentage = rep()
pathway_name1 = names(Matrix_gene_sets_morethan2)
for (i in 1:length(Matrix_gene_sets_morethan2)){
   pathway_lengths2 = c(pathway_lengths2, length(Matrix_gene_sets_morethan2[[i]]))
   percentage = c(percentage, 
                  length(Matrix_gene_sets_morethan2[[i]])/
                    length(Matrix_gene_sets_unconverted_morethan2[[i]]))
}

#save converted c2 gene files
saveRDS(Matrix_gene_sets_converted, 'Matrix_c2gene_sets_converted.rds')
#save converted c2 gene files taking out the empty ones and the ones less than 2 genes
saveRDS(Matrix_gene_sets_morethan2, 'Matrix_gene_sets_morethan2.rds')
###write list to csv
Matrix_gene_sets_morethan2.df <- do.call("rbind", sapply(Matrix_gene_sets_morethan2, as.data.frame)) 
write.csv(Matrix_gene_sets_morethan2.df, file = "Matrix_gene_sets_morethan2.csv")

#correct dataframe containing the pathway names and number of proteins each
pathway_name1 = names(Matrix_gene_sets_morethan2)
pathways_df1 = data.frame(pathway = pathway_name1, protein_num = pathway_lengths2, coverage = percentage)
write_csv(pathways_df1, 'pathway_df.csv')

#maximum coverage
max(pathways_df1['coverage'])

######### only look at pathways containing more than 5 genes
pathway_to_delete = rep()
pathway_lengths1 = rep()
pathway_lengths2 = rep()
for (i in 1:length(Matrix_gene_sets_non_null)){ 
  if (length(Matrix_gene_sets_non_null[[i]]) < 6) {
    pathway_to_delete = c(pathway_to_delete, i)
  }
  pathway_lengths1 = c(pathway_lengths1, length(Matrix_gene_sets_non_null[[i]]))
}
Matrix_gene_sets_morethan5 = Matrix_gene_sets_non_null[-pathway_to_delete]
#get the unconverted pathway matrix to be of the same length as Matrix_gene_sets_morethan2
Matrix_gene_sets_c2_nonnull5 = Matrix_gene_sets_c2[-empty_pathways]
Matrix_gene_sets_unconverted_morethan5 = Matrix_gene_sets_c2_nonnull5[-pathway_to_delete]

#calculate the percentage of genes in each pathway
percentage = rep()
pathway_name1 = rownames(Matrix_gene_sets_morethan5)
for (i in 1:length(Matrix_gene_sets_morethan5)){
  pathway_lengths2 = c(pathway_lengths2, length(Matrix_gene_sets_morethan5[[i]]))
  percentage = c(percentage, 
                 length(Matrix_gene_sets_morethan5[[i]])/
                   length(Matrix_gene_sets_unconverted_morethan5[[i]]))
}
#save converted c2 gene files
saveRDS(Matrix_gene_sets_converted, 'Matrix_c2gene_sets_converted.rds')
#save converted c2 gene files taking out the empty ones and the ones less than 2 genes
saveRDS(Matrix_gene_sets_morethan5, 'Matrix_gene_sets_morethan5.rds')
###write list to csv
Matrix_gene_sets_morethan5.df <- do.call("rbind", sapply(Matrix_gene_sets_morethan2, as.data.frame)) 
write.csv(Matrix_gene_sets_morethan5.df, file = "Matrix_gene_sets_morethan5.csv")

#correct dataframe containing the pathway names and number of proteins each
pathway_name1 = names(Matrix_gene_sets_morethan5)
pathways_df5 = data.frame(pathway = pathway_name1, protein_num = pathway_lengths2, coverage = percentage)
write_csv(pathways_df5, 'pathway_df5.csv')

#maximum coverage
max(pathways_df5['coverage'])


############
#converting mice gene to human protein in the c5's pathways
unconverted = 0
converted = 0
empty_pathways = rep()
pathway_lengths = rep()
Matrix_gene_sets_c5 = readRDS('c5_gene_sets_beforeconvert.rds')
Matrix_gene_sets_converted = vector(mode = "list", length = length(Matrix_gene_sets_c5))
for (i in 1:length(Matrix_gene_sets_c5)){
  for (micegene in Matrix_gene_sets_c5[[i]]){
    if (micegene != ""){
      if (micegene %in% toupper(full_converter$mice_gene)){
        human = full_converter$human_gene_clean[toupper(full_converter$mice_gene) == micegene]
        Matrix_gene_sets_converted[[i]] = c(Matrix_gene_sets_converted[[i]], human)
        converted = converted + 1
      } else {
        unconverted = unconverted + 1
      }
    }
  }
  #get a list of pathway index if the pathway contains two or less genes
  #if (is.null(Matrix_gene_sets_converted[[i]])) {empty_pathways = c(empty_pathways, i)}
  if (length(Matrix_gene_sets_converted[[i]]) < 3) {
    empty_pathways = c(empty_pathways, i)
  }
  pathway_lengths = c(pathway_lengths, length(Matrix_gene_sets_converted[[i]]))
}

#a dataframe containing the pathway names and number of proteins each
# pathways_df = data.frame(pathyway = pathway_name, protein_num = pathway_lengths)
# write_csv(pathways_df, 'pathway_df.csv')

names(Matrix_gene_sets_converted) <- names(Matrix_gene_sets_c5)

#Taking out the empty ones and the ones less than 2 genes
Matrix_gene_sets_non_null = Matrix_gene_sets_converted[-empty_pathways] #delete the ones with empty pathways
pathway_to_delete = rep()
pathway_lengths1 = rep()
pathway_lengths2 = rep()
for (i in 1:length(Matrix_gene_sets_non_null)){ 
  if (length(Matrix_gene_sets_non_null[[i]]) < 2) {
    pathway_to_delete = c(pathway_to_delete, i)
  }
  pathway_lengths1 = c(pathway_lengths1, length(Matrix_gene_sets_non_null[[i]]))
}
Matrix_gene_sets_morethan2 = Matrix_gene_sets_non_null#[-pathway_to_delete]

#get the unconverted pathway matrix to be of the same length as Matrix_gene_sets_morethan2
Matrix_gene_sets_c5_nonnull = Matrix_gene_sets_c5[-empty_pathways]
Matrix_gene_sets_unconverted_morethan2 = Matrix_gene_sets_c5_nonnull#[-pathway_to_delete]
#write to csv Matrix_gene_sets_c5_nonnull.df
write.table(Matrix_gene_sets_c5_nonnull, file = "Matrix_gene_sets_c5_nonnull.csv")
Matrix_gene_sets_c5_nonnull.df <- do.call("rbind", sapply(Matrix_gene_sets_c5_nonnull, as.data.frame)) 
write.csv(Matrix_gene_sets_c5_nonnull.df, file = "Matrix_gene_sets_c5_nonnull.csv")
pathway_name_c5_nonnull = names(Matrix_gene_sets_c5_nonnull)
write.csv(pathway_name_c5_nonnull, "pathway_name_c5_nonnull.csv")
#calculate the percentage of genes successfully converted to human in each pathway:
#num of proteins converted found in the pathway/num of protein in the pathway total before conversion
percentage = rep()
pathway_name1 = names(Matrix_gene_sets_morethan2)
for (i in 1:length(Matrix_gene_sets_morethan2)){
  pathway_lengths2 = c(pathway_lengths2, length(Matrix_gene_sets_morethan2[[i]]))
  percentage = c(percentage, 
                 length(Matrix_gene_sets_morethan2[[i]])/
                   length(Matrix_gene_sets_unconverted_morethan2[[i]]))
}

#save converted c5 gene files
saveRDS(Matrix_gene_sets_converted, 'Matrix_c5gene_sets_converted.rds')
#save converted c5 gene files taking out the empty ones and the ones less than 2 genes
saveRDS(Matrix_gene_sets_morethan2, 'Matrix_gene_sets_morethan2_c5.rds')
###write list to csv
Matrix_gene_sets_morethan2.df <- do.call("rbind", sapply(Matrix_gene_sets_morethan2, as.data.frame)) 
write.csv(Matrix_gene_sets_morethan2.df, file = "Matrix_gene_sets_morethan2_c5.csv")

#correct dataframe containing the pathway names and number of proteins each
pathway_name1 = names(Matrix_gene_sets_morethan2)
pathways_df1 = data.frame(pathway = pathway_name1, protein_num = pathway_lengths2, coverage = percentage)
write_csv(pathways_df1, 'pathway_df_c5.csv')

#maximum percentage of genes converted to human gene names
max(pathways_df1['coverage'])


