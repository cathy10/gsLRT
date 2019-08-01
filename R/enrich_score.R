#'
#'
#'
#'
enrich_score <- function(stat_mat, gene_sets, trim = T, min_size = 5, stat = "bb"){

  require(magrittr)
  require(plyr)

  gene_names <- colnames(stat_mat)

  if(trim){
    gene_sets %<>% plyr::llply(function(gs){
      ind <- gs %in% gene_names
      if(sum(ind) < min_size){return(NULL)}else{return(gs[ind])}
    })
    gene_sets[sapply(gene_sets, is.null)] <- NULL
  }

  if(stat == "robust_chisq"){
    gs_sizes <- plyr::laply(gene_sets, function(gs){length(gs)})
    emp_dist <- plyr::alply(1:max(gs_sizes), 1, function(s){
      plyr::aaply(1:10000, 1, function(i){
        median(rgamma(s, 1/2, 1/2))
      })
    })
  }

  plyr::laply(gene_sets, function(gs){

    gs_ind <- ifelse(gene_names %in% gs, 1, 0)
    gs_bool <- gs_ind == 1
    card <- length(gs)

    if(stat == "bb"){

      plyr::aaply(1:nrow(stat_mat), 1, function(i){

        c_i <- stat_mat[i,]
        ord <- order(c_i, decreasing = T)
        c_ri <- c_i[ord]
        gs_rind <- gs_ind[ord]
        rw <- (cumsum(abs(c_ri)*gs_rind) / sum(abs(c_ri)*gs_rind)) %>%
          magrittr::subtract(cumsum(1-gs_rind)/(length(gene_names) - card))
        if(any(is.na(rw))){
          print("found na")
          return(0)
        } else {
          rmax <- max(rw, na.rm = T)
          rmin <- min(rw, na.rm = T)
          if(abs(rmax) > abs(rmin)){
            return(rmax)
          } else {
            return(rmin)
          }
        }
      })

    } else if(stat == "chisq"){

      sum(stat_mat[, gs_bool]) %>%
          pgamma(shape = card / 2, rate = 1 / 2, log.p = T) %>%
          qgamma(shape = 1 / 2, rate = 1 / 2, log.p = T)

    } else if(stat == "hg"){

      qchisq <- qgamma(0.05, 0.5, 0.5)
      yy <- sum(stat_mat[, gs_bool] > qchisq)
      yn <- sum(stat_mat[, !gs_bool] > qchisq)
      ny <- sum(stat_mat[, gs_bool] < qchisq)
      nn <- sum(stat_mat[, !gs_bool] < qchisq)
      matrix(c(yy, yn, ny, nn), nrow = 2, ncol = 2) %>%
        fisher.test(alternative = "greater") %>%
        magrittr::extract2("p.value")

    }

  }) %>%
    t() %>%
    magrittr::set_colnames(names(gene_sets))

}

