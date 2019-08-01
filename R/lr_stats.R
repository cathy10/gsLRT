#'
#'
#'
#'
lr_stats <- function(X, y, z,
                     interaction_term = NULL, stat = "wald"){
  require(magrittr)
  cnames <- colnames(X)
  X <- cbind(rep(1, length(y)), X) %>% magrittr::set_colnames(c("Intercept", cnames))
  X %<>% apply(2, function(x){
    if(length(unique(x)) > 2){
      (x - mean(x)) / sd(x)
    } else {
      x
    }
  })
  if(is.null(interaction_term)){
    df_j <- data.frame(g = (y - mean(y, na.rm=T)) / sd(y, na.rm=T),
                       X,
                       class = z)
    formula_j0 <- paste0(paste(colnames(X),
                                collapse = " + ")) %>%
                 paste0("class ~ 0 + ", .)
    mod_j0 <- glm(formula_j0, family = "binomial", data = df_j)
    fv_j0 <- fitted.values(mod_j0)
    formula_j1 <- paste0(paste(colnames(X),
                                collapse = " + ")) %>%
                 paste0("class ~ 0 + g + ", .)
    mod_j1 <- glm(formula_j1, family = "binomial", data = df_j)
    fv_j1 <- fitted.values(mod_j1)
    std_err <- summary(mod_j1)[["coefficients"]][, "Std. Error"]["g"]
  } else {
    df_j0 <- data.frame(g = (y - mean(y, na.rm=T)) / sd(y, na.rm=T),
                       X,
                       class = z)
    formula_j0 <- paste0(paste(colnames(X),
                                collapse = " + ")) %>%
                 paste0("class ~ 0 + g + ", .)
    mod_j0 <- glm(formula_j0, family = "binomial", data = df_j0)
    fv_j0 <- fitted.values(mod_j0)
    df_j1 <- df_j0 %>%
              dplyr::mutate(interaction = df_j0[["g"]]*df_j0[[interaction_term]])
    formula_j1 <- paste0(paste(colnames(X),
                                collapse = " + ")) %>%
                 paste0("class ~ 0 + g + interaction + ", .)
    mod_j1 <- glm(formula_j1, family = "binomial", data = df_j1)
    fv_j1 <- fitted.values(mod_j1)
    std_err <- summary(mod_j1)[["coefficients"]][, "Std. Error"]["interaction"]
  }
  if(stat == "wald"){
    c_j <- (log(fv_j1 / (1-fv_j1)) - log(fv_j0 / (1-fv_j0))) / std_err
  } else if(stat == "bartlett"){
    c_j <- 2*(z*log(fv_j1) + (1 - z)*log(1 - fv_j1)) -
            2*(z*log(fv_j0) + (1 - z)*log(1 - fv_j0))
    alpha_n <- (nrow(X)/2)*(sum(hatvalues(mod_j1)^2) - sum(hatvalues(mod_j0)^2))
    c_j <- c_j / (1 + (alpha_n / nrow(X)))
  }
  return(c_j)
}
