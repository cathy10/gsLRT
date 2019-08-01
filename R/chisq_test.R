#'
#'
#'
chisq_test <- function(X, gs_list){
  require(plyr)
  require(magrittr)

  plyr::adply(colnames(X), 1, function(gs){

    stat <- sum(X[, gs])

    data.frame(GS = gs,
               stat = stat,
               p = 1 - exp(pgamma(stat, nrow(X)/2, 1/2, log.p = T)),
               genes = gs_list[[gs]] %>%
                 sort() %>%
                 paste(collapse = ", "),
               size = length(gs_list[[gs]]))

  }) %>%
    dplyr::mutate(q = p.adjust(p, method = "BH")) %>%
    plyr::arrange(p)

}
