plot_heatmap <- function(X, col_annotation = NULL, row_annotation = NULL,
                         annotation_colors = NULL,
                         arrange_by = NULL, col_range = c(-1, 1)){

  if(!is.null(col_annotation)){
    c_annot <- do.call(cbind, col_annotation) %>%
                as.data.frame() %>%
                magrittr::set_colnames(names(col_annotation))  %>%
                dplyr::mutate(row_name = row.names(X))
    if(!is.null(arrange_by)){
      c_annot %<>% dplyr::arrange_at(which(colnames(.) %in% arrange_by))
      X <- X[, c_annot[["row_name"]]]
    }
    c_annot %<>%
      magrittr::set_rownames(magrittr::extract2(., "row_name")) %>%
      dplyr::select(-row_name)

  } else {
    c_annot <- NULL
  }

  if(!is.null(row_annotation)){
    r_annot <- do.call(cbind, row_annotation) %>%
      as.data.frame() %>%
      magrittr::set_colnames(names(row_annotation)) %>%
      dplyr::mutate(row_name = row.names(X))
    if(!is.null(arrange_by)){
      r_annot %<>% dplyr::arrange_at(which(colnames(.) %in% arrange_by))
      X <- X[r_annot[["row_name"]], ]
    }
    r_annot %<>%
      magrittr::set_rownames(magrittr::extract2(., "row_name")) %>%
      dplyr::select(-row_name)
  } else {
    r_annot <- NULL
  }


  cluster_cols <- ifelse(any(arrange_by %in% names(col_annotation)), F, T)
  cluster_rows <- ifelse(any(arrange_by %in% names(row_annotation)), F, T)

  pheatmap::pheatmap(X,
                     annotation_col = c_annot,
                     annotation_row = r_annot,
                     breaks=seq(col_range[1],col_range[2],length.out=101),
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, name = "RdBu")))(100),
                     show_rownames = F, show_colnames = F,
                     annotation_colors = annotation_colors,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols)

}
