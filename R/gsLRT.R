#'
#'
#'
#'
gsLRT <- function(y, X, Z, gene_sets, interaction_term = NULL, level = "gene", n_perm = 100, do_par = T){
  #
  require(magrittr)
  require(plyr)
  require(doParallel)
  doParallel::registerDoParallel(cores = 4)
  #
  perm_results <-
    plyr::aaply(1:(n_perm + 1), 1, function(i){
      #
      if(i == 1){
        r_ind <- 1:nrow(Z)
      } else {
        r_ind <- sample(1:nrow(Z), nrow(Z))
      }
      #
      cnames <- colnames(X)
      #
      X <- cbind(rep(1, nrow(Z)), X) %>% magrittr::set_colnames(c("Intercept", cnames))
      X %<>% apply(2, function(x){
        if(length(unique(x)) > 2){
          (x - mean(x)) / sd(x)
        } else {
          x
        }
      })
      #
      plyr::aaply(1:ncol(Z), 1, function(j){
        #
        if(is.null(interaction_term)){
          #
          z <- Z[r_ind, j]
          #
          df_j <- data.frame(g = (z - mean(z, na.rm=T)) / sd(z, na.rm=T),
                             X,
                             class = y)
          formula_j0 <- paste0(paste(colnames(X),
                                     collapse = " + ")) %>%
            paste0("class ~ 0 + ", .)
          mod_j0 <- glm(formula_j0, family = "binomial", data = df_j)
          #mod_j0 <- RcppNumerical::fastLR(X, y)
          fv_j0 <- mod_j0[["fitted.values"]]
          formula_j1 <- paste0(paste(colnames(X),
                                     collapse = " + ")) %>%
            paste0("class ~ 0 + g + ", .)
          mod_j1 <- glm(formula_j1, family = "binomial", data = df_j)
          #mod_j1 <- RcppNumerical::fastLR(cbind(X, df_j[["g"]]), y)
          fv_j1 <- mod_j1[["fitted.values"]]
          std_err <- summary(mod_j1)[["coefficients"]][, "Std. Error"]["g"]
        } else {
          #
          z <- Z[, j]
          #
          df_j0 <- data.frame(g = (z - mean(z, na.rm=T)) / sd(z, na.rm=T),
                              X,
                              class = y)
          formula_j0 <- paste0(paste(colnames(X),
                                     collapse = " + ")) %>%
            paste0("class ~ 0 + g + ", .)
          mod_j0 <- glm(formula_j0, family = "binomial", data = df_j0)
          fv_j0 <- fitted.values(mod_j0)
          df_j1 <- df_j0 %>%
            dplyr::mutate(interaction = (df_j0[["g"]]*df_j0[[interaction_term]])[r_ind])
          formula_j1 <- paste0(paste(colnames(X),
                                     collapse = " + ")) %>%
            paste0("class ~ 0 + g + interaction + ", .)
          mod_j1 <- glm(formula_j1, family = "binomial", data = df_j1)
          fv_j1 <- fitted.values(mod_j1)
        }
        #
        if(level == "gene"){
          c_j <- 2*(y*log(fv_j1) + (1 - y)*log(1 - fv_j1)) -
            2*(y*log(fv_j0) + (1 - y)*log(1 - fv_j0))
          alpha_n <- (nrow(X)/2)*(sum(hatvalues(mod_j1)^2) - sum(hatvalues(mod_j0)^2))
          c_j <- c_j / (1 + (alpha_n / nrow(X)))
          return(sum(c_j))
        } else if(level == "sample"){
          c_j <- (log(fv_j1 / (1-fv_j1)) - log(fv_j0 / (1-fv_j0))) / std_err
          return(c_j)
        }
      })

    }, .parallel = do_par, .progress = "text") %>%
    magrittr::set_colnames(colnames(Z))
  #
  if(level == "gene"){
    #
    perm_mean <- apply(perm_results[2:nrow(perm_results), ], 1, mean)
    perm_sd <- apply(perm_results[2:nrow(perm_results), ], 1, sd)
    rand_mean <- mean(perm_results[1, ])
    rand_sd <- sd(perm_results[1, ])
    #
    restandardized_results <-
      plyr::adply(names(gene_sets), 1, function(ngs){

        gs <- gene_sets[[ngs]]
        if(length(gs) > 1){
          S <- mean(perm_results[1, gs])
          S_star <- rowMeans(perm_results[2:nrow(perm_results), gs])
          S_star_star <- (S_star - perm_mean) / perm_sd
          data.frame(GS = ngs,
                     genes = gs %>%
                       sort() %>%
                       paste(collapse = ", "),
                     size = length(gs),
                     s = S,
                     p = (sum(S_star > S) + 1) / (length(S_star_star) + 1))
        }

      }) %>%
      plyr::arrange(p) %>%
      dplyr::mutate(q = p)
    return(restandardized_results)
  } else if(level == "sample") {
    #
    perm_mean <- apply(perm_results[2:nrow(perm_results), ,], c(1,3), mean)
    perm_sd <- apply(perm_results[2:nrow(perm_results), ,], c(1,3), sd)
    rand_mean <- apply(perm_results[1, ,], 2, mean)
    rand_sd <- apply(perm_results[1, ,], 2, sd)
    #
    restandardized_results <-
      plyr::adply(names(gene_sets), 1, function(ngs){

        gs <- gene_sets[[ngs]]
        if(length(gs) > 1){
          plyr::adply(1:dim(perm_results)[3], 1, function(i){
            S <- mean(perm_results[1, gs, i])
            S_star <- rowMeans(perm_results[2:nrow(perm_results), gs, i])
            S_star_star <- (S_star - perm_mean[, i]) / perm_sd[, i]
            data.frame(sample_id = i,
                       GS = ngs,
                       genes = gs %>%
                         sort() %>%
                         paste(collapse = ", "),
                       size = length(gs),
                       s = S,
                       p = (sum(S_star > S) + 1) / (length(S_star_star) + 1))

          })
        }

      }) %>%
      plyr::arrange(p) %>%
      dplyr::mutate(q = p)
    return(restandardized_results)
  }
}
