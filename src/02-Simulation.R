# Simulate data according to linear mixed effects model
## Load libraries
n_sims <- 25
sim_df <- plyr::adply(1:n_sims, 1, function(sim){
  library(plyr)
  library(tidyverse)
  library(magrittr)
  library(MASS)
  library(lme4)
  ## Set params
  n <- 152 # number of samples
  p <- 142 # number of genes
  d <- 3 # number of continuous covariates
  b <- 3 # number of binary covariates
  m <- 20 # number of gene sets
  tau_sq <- 1 # signal variance
  ## Create fake genes and gene sets
  gene_names <- sample(LETTERS, 6*p, replace = T) %>%
    matrix(ncol=6) %>%
    apply(1, paste, collapse="")
  #
  gene_sets <- split(gene_names, sample(1:m, length(gene_names), 
                                        replace=T, 
                                        prob = runif(m) %>%
                                          magrittr::divide_by(., sum(.))))
  gene_sets[sapply(gene_sets, function(s){length(s) < 4})] <- NULL
  ## Set values for design matrix
  X <- matrix(rnorm(n*d, sd=1), nrow = n, ncol = d) %>% 
    magrittr::set_colnames(paste0("c_var", 1:d)) %>%
    cbind(matrix(sample(c(0, 1), n*b, replace = T), nrow = n, ncol = b) %>%
            magrittr::set_colnames(paste0("b_var", 1:b)))
  betas <- rnorm(b + d)
  ## Set gene sets that should have enrichment
  gs_sizes <- plyr::laply(gene_sets, function(x){
    length(x)
  })
  e_gene_sets <- sample(names(gene_sets)[gs_sizes > 5], 5)
  ## Create gene expression matrix
  e_genes <- plyr::alply(e_gene_sets, 1, function(x){
    gene_sets[[x]]
  }) %>%
    unlist()
  #
  G <- plyr::llply(gene_sets, function(gs){
    df <- length(gs) + length(gs)
    cor_strength <- runif(1)
    S0_inv <- rep(cor_strength, length(gs)^2) %>%
      matrix(nrow = length(gs)) %>%
      magrittr::add(diag(1 - cor_strength, nrow = length(gs), ncol = length(gs))) %>%
      MASS::ginv()
    Z <- t(chol(S0_inv))%*%matrix(rnorm(length(gs)*df), nrow=length(gs)) ## Invert matrix
    Sigma_inv <- Z%*%t(Z)
    Sigma <- MASS::ginv(Sigma_inv)
    L <- chol(Sigma)
    #L <- chol(diag(rep(1, length(gs))))
    (t(L)%*%matrix(rnorm(length(gs)*n), nrow = length(gs))) %>%
      t() %>%
      magrittr::set_colnames(gs)
    
  }) %>%
    do.call(cbind, .) %>%
    apply(2, function(x){(x - mean(x))/sd(x)})
  #
  gammas <- plyr::adply(e_gene_sets, 1, function(gs_i){
    
    genes <- gene_sets[[gs_i]]
    data.frame(gene_name = genes,
               coef = rnorm(length(genes), sd = sqrt(tau_sq)))
    
  }) %>%
    dplyr::mutate(ph = 1) %>%
    reshape2::acast(ph~gene_name, value.var = 'coef', mean)
  
  Z <- G[, e_genes] %*%
    gammas[, e_genes] %>%
    magrittr::add(X%*%betas) %>%
    exp() %>%
    magrittr::divide_by(., 1 + .) %>%
    rbinom(n = length(.), size = 1, prob = .) %>%
    as.integer()
  
  ## Run MMA
  presults <- gsHyp::gsHyp(Z, X, G, gene_sets = gene_sets, n_perm = 200, do_par = T)
  #
  romer_df <- limma::romer(t(G), 
                           gene_sets %>%
                             plyr::llply(function(gs){
                               which(colnames(G) %in% gs)
                             }), 
                           design = cbind(Intercept = 1, X, Group = Z),
                           nrot = 200) %>%
    data.frame() %>%
    dplyr::mutate(GS = row.names(.), p = Mixed)
  #
  gs_hyp_fp <- c()
  gs_hyp_tp <- c()
  romer_fp <- c()
  romer_tp <- c()
  #
  pos <- sum(e_gene_sets %in% presults$GS)
  neg <- length(unique(presults$GS)) - pos
  #
  for(alpha in seq(log(0.005), log(1), length.out = 50) %>% exp()){
    gs_hyp_fp %<>% c(sum(!(presults$GS[presults$p < alpha] %in% e_gene_sets)) / neg)
    gs_hyp_tp %<>% c(sum((e_gene_sets %in% presults$GS[presults$p < alpha])) / pos)
    romer_fp %<>% c(sum(!(romer_df$GS[romer_df$p < alpha] %in% e_gene_sets)) / neg)
    romer_tp %<>% c(sum((e_gene_sets %in% romer_df$GS[romer_df$p < alpha])) / pos)
  }
  #
  data.frame(sim = 1,
             alpha = seq(log(0.005), log(1), length.out = 50) %>% exp(),
             gs_hyp_fp = gs_hyp_fp,
             gs_hyp_tp = gs_hyp_tp,
             romer_fp = romer_fp,
             romer_tp = romer_tp)
  
}, .progress = "text")
#
sim_df %>%
  readr::write_tsv("results/tables/simulation_norm.txt")
#
rbind(sim_df %>% 
        dplyr::select(sim, alpha, gs_hyp_tp, gs_hyp_fp) %>% 
        dplyr::rename(tp = gs_hyp_tp, fp = gs_hyp_fp) %>%
        dplyr::mutate(type = "gsHyp"),
      sim_df %>% 
        dplyr::select(sim, alpha, romer_tp, romer_fp) %>% 
        dplyr::rename(tp = romer_tp, fp = romer_fp) %>%
        dplyr::mutate(type = "ROAST")) %>%
  dplyr::group_by(alpha, type) %>%
  dplyr::summarise(tp = mean(tp), fp = mean(fp)) %>%
  ggplot2::ggplot() +
  geom_point(aes(x = alpha, y = tp, colour = type, group = type)) +
  geom_line(aes(x = alpha, y = tp, colour = type, group = type)) +
  ggthemes::theme_few() +
  ggthemes::scale_color_calc(name = "Method") +
  labs(x = "p-value threshold", y = "True positive rate")
ggplot2::ggsave("results/figures/compare_tpr.jpeg", height = 5, width = 6)
#
rbind(sim_df %>% 
        dplyr::select(sim, alpha, gs_hyp_tp, gs_hyp_fp) %>% 
        dplyr::rename(tp = gs_hyp_tp, fp = gs_hyp_fp) %>%
        dplyr::mutate(type = "gsHyp"),
      sim_df %>% 
        dplyr::select(sim, alpha, romer_tp, romer_fp) %>% 
        dplyr::rename(tp = romer_tp, fp = romer_fp) %>%
        dplyr::mutate(type = "ROAST")) %>%
  dplyr::group_by(alpha, type) %>%
  dplyr::summarise(tp = mean(tp), fp = mean(fp)) %>%
  ggplot2::ggplot() +
  geom_point(aes(x = alpha, y = fp, colour = type, group = type)) +
  geom_line(aes(x = alpha, y = fp, colour = type, group = type)) +
  ggthemes::theme_few() +
  ggthemes::scale_color_calc(name = "Method") +
  labs(x = "p-value threshold", y = "False positive rate")
ggplot2::ggsave("results/figures/compare_fpr.jpeg", height = 5, width = 6)
#
rbind(sim_df %>% 
        dplyr::select(sim, alpha, gs_hyp_tp, gs_hyp_fp) %>% 
        dplyr::rename(tp = gs_hyp_tp, fp = gs_hyp_fp) %>%
        dplyr::mutate(type = "gsHyp"),
      sim_df %>% 
        dplyr::select(sim, alpha, romer_tp, romer_fp) %>% 
        dplyr::rename(tp = romer_tp, fp = romer_fp) %>%
        dplyr::mutate(type = "ROAST")) %>%
  dplyr::group_by(alpha, type) %>%
  dplyr::summarise(tp = mean(tp), fp = mean(fp)) %>%
  ggplot2::ggplot() +
  geom_point(aes(x = fp, y = tp, colour = type, group = type)) +
  geom_line(aes(x = fp, y = tp, colour = type, group = type)) +
  ggthemes::theme_few() +
  ggthemes::scale_color_calc(name = "Method") +
  labs(x = "False positive rate", y = "True positive rate")
ggplot2::ggsave("results/figures/compare_auc.jpeg", height = 5, width = 6)
  


