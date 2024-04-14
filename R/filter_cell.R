#' Filter cells
#' @description
#' Filter the cells that are non-significant perturbed.
#'
#' @param count A seurat object that has been normalized.
#' @param sgrna A cell-perturbation data frame.
#' @param filter Ther ratio of filtering cells that are perturbed non-significant. Default: 0.5.
#' @param npcs A number of PCs in PCA. Default: 20 (Same as number of factor)
#' @param nfeatures A parameter for extracting high-variable genes. Default: 3000.
#' @param ctrl The name of control perturbation. Default: NT
#'
#' @return A filtered seurat object.
#' @export
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' rna_filter<-filter_cell(rna,sgrna,filter=0.5,ctrl='NT')
filter_cell <- function(count, sgrna, filter=0.5, npcs=20, nfeatures=3000, ctrl='NT'){

  normal_data<-LayerData(count, assay = "RNA", layer = 'data')

  pert <- as.character(unique(sgrna$perturbation))
  pert1 <- pert[-which(pert == ctrl)]

  if(sum(pert1 %in% rownames(normal_data)) == 0){
    stop('The perturbations do not in the rownames of count data. Please check and convert gene type into symbol!')
  }

  perturbation_data <- sgrna
  for(i in 1:length(pert1)){
    a<-(normal_data[,which(colnames(normal_data) %in% sgrna$cell[which(sgrna$perturbation == pert1[i])])])
    a<-a[which(rownames(a) == pert1[i]),]
    a<-a[order(a,decreasing = T)]
    a<-a[1:(round(length(a) * filter))]
    perturbation_data<-perturbation_data[-which(perturbation_data$cell %in% names(a)),]
  }
  perturbation_data<-perturbation_data[-which(perturbation_data$perturbation==ctrl),]

  combined_obj_nt<-count[,which(colnames(count) %in% sgrna$cell[which(sgrna$perturbation == ctrl)])]
  combined_obj_nt <- NormalizeData(combined_obj_nt) %>% FindVariableFeatures(nfeatures = nfeatures) %>% ScaleData()
  combined_obj_nt<-RunPCA(combined_obj_nt,npcs = npcs)
  combined_obj_nt<-FindNeighbors(combined_obj_nt, dims = 1:npcs) %>% FindClusters()
  nt_cell<-colnames(combined_obj_nt[,which(combined_obj_nt@active.ident == '0')])

  rna_filter<-count[,which(colnames(count) %in% c(perturbation_data$cell, nt_cell))]#7395cells
  rna_filter<-NormalizeData(rna_filter) %>% FindVariableFeatures(nfeatures = nfeatures) %>% ScaleData(do.center = F) #注意do.center = F

  return(rna_filter)
}
