#' Calculate perturbmeta score for each cluster
#' @description
#' Euclidean distance is used to calculate the score of the perturbation to each factor, and then
#' Monte Carlo permutation is used to test the P-value, and multiple tests are used to correct the P-value.
#'
#' @param nmf The result from function perturbmeta.
#' @param sgrna Ther cell-perturbation data frame.
#' @param cluster The number of clusters that used for perturbmeta.
#' @param factor The number of factor in NMF. Default: 20.
#' @param nperms The number of permutations. Default: 1000.
#'
#' @return A list contain perturbmeta_score and p value.
#' @export
#' @import dplyr
#' @import tibble
#'
#' @examples
#' score_res<-perturbmeta_score(res_nmf,sgrna,cluster=5,factor=20,nperms=1000)
perturbmeta_score <- function(nmf, sgrna, cluster=5, factor=20, nperms=1000){
  PerturbMeta_score_list <- list()
  for (i in 1:cluster){
    expr = t(nmf[[i]]@fit@H) %>% data.frame()
    metadata2 = sgrna[sgrna$cell %in% rownames(expr),]
    expr$perturbation <- metadata2$perturbation[match(rownames(expr),metadata2$cell)]
    PerturbMeta_score <- expr %>% group_by(perturbation) %>% summarise(across(where(is.numeric), ~mean((.^2)))) %>% column_to_rownames("perturbation") %>% scale()
    colnames(PerturbMeta_score) <- paste0("factor",seq(1:factor))
    PerturbMeta_score_list[[i]] <- PerturbMeta_score
  }

  results <- list()
  n_perms <- nperms
  for (i in 1:cluster){
    cluster_cell_embedding <- nmf[[i]]@fit@H
    rownames(cluster_cell_embedding) <- paste0("factor",1:factor)
    group <- sgrna$perturbation[match(colnames(cluster_cell_embedding),sgrna$cell)]
    factor_results <- list()
    for (j in rownames(cluster_cell_embedding)){

      cell_set_score <- cluster_cell_embedding[j, ]
      perturbation_results <- list()
      for (perturbation in unique(group)){

        n_perms <- nperms
        observed_score <- mean(cell_set_score[group == perturbation]) - mean(cell_set_score[group != perturbation])
        perm_scores <- numeric(nperms)

        set.seed(123)
        for (n in 1:nperms) {

          perm_group <- sample(group)

          perm_score <- mean(cell_set_score[perm_group == perturbation]) - mean(cell_set_score[perm_group != perturbation])

          perm_scores[n] <- perm_score
        }

        perm_mean <- mean(perm_scores)
        perm_sd <- sd(perm_scores)

        z_score <- (observed_score - perm_mean) / perm_sd

        p_value <- 2 * (1 - pnorm(abs(z_score)))
        perturbation_results[[perturbation]] <- p_value
      }
      factor_results[[j]] <- p.adjust(unlist(perturbation_results),method = "bonferroni")
    }
    results[[i]] <- factor_results
  }
  results_df_bonferroni <- lapply(results, function(x) {do.call("rbind",x)})

  return(list(PerturbMeta_score_list, results_df_bonferroni))
}
