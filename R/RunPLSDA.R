#' Run PLS-DA
#' @description
#' Perform PLS-DA for each cluster, then extract high-contribution genes for further analysis.
#'
#' @param count A count data from selected cluster.
#' @param sgrna A cell-perturbation data frame.
#' @param VIP The threshold of extracting features from PLS-DA. Default: 1
#' @param ncomp A number of component in PLS-DA. Default: 20 (Same as number of factor)
#'
#' @return The list contain cell-embedding, feature-loading and high-contribution genes.
#' @export
#' @import mixOmics
#' @import RVAideMemoire
#'
#' @examples
#' plsda_res<-RunPLSDA(count,sgrna,ncomp=20,VIP=1)
RunPLSDA <- function(count, sgrna, VIP=1, ncomp=20){

  out<-count
  out<-t(out)

  group<-sgrna
  pert<-as.character(unique(group$perturbation))
  group$perturbation<-as.character(group$perturbation)

  group1<-group
  group1<-group1[which(group1$cell %in% rownames(out)),]
  df_plsda<-plsda(out, group1$perturbation,ncomp = ncomp)
  embed<-data.frame(df_plsda$variates$X)
  load<-data.frame(df_plsda$loadings$X)

  da2<-PLSDA.VIP(df_plsda, graph = F)
  da2<-data.frame(da2$tab)
  da2<-cbind(rownames(da2),da2)
  da2<-da2[da2$VIP>VIP,]

  return(list(embed, load, da2))
}
