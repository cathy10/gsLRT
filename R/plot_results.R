#'
#'
#'
#'
plot_results <- function(lm_res, stats, gs_list){

  top_gs <- lm_res %>%
              dplyr::top_n(8, 1/p)
  sum_df <- stats %>%
              reshape2::melt(varnames = c("sample", "gene"),
                             value.name = "stat") %>%
              dplyr::group_by(., gene) %>%
              dplyr::summarise(sum = sum(stat)) %>%
              dplyr::ungroup()

  p_df <- plyr::adply(1:nrow(top_gs), 1, function(i){

    data.frame(GS = top_gs[["GS"]][i],
               q = top_gs[["q"]][i],
               p = top_gs[["p"]][i],
               gene = gs_list[[top_gs[["GS"]][i]]])

  }) %>%
    dplyr::left_join(sum_df) %>%
    dplyr::mutate(GS = factor(GS, levels = unique(GS[order(p)]))) %>%
    dplyr::group_by(GS) %>%
    dplyr::top_n(10, sum) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(GS, sum) %>%
    dplyr::mutate(order = row_number())

  p_df %>%
    ggplot2::ggplot() +
    geom_bar(aes(x = order, y = sum, fill = GS), stat = "identity") +
    facet_grid(GS~., scales = "free", switch = "both") +
    coord_flip() +
    scale_x_continuous(
      breaks = p_df$order,
      labels = p_df$gene,
      expand= c(0, 0),
      position = "top"
    ) +
    ggthemes::theme_few() +
    ggthemes::scale_fill_colorblind() +
    guides(fill = F) +
    theme(strip.text.y = element_text(angle = 180, size = 6),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 5, angle = 30, hjust = 1)) +
    labs(y = "LLR")

}
