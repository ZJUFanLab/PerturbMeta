#' Calculate PLS-DA and NMF for each cluster
#' @description
#' Filter the cells that are non-significant perturbed, then carry out NMF, extract the cell-embedding and feature-loading for
#' re-clustering analysis. And perform PLS-DA for each cluster, then extract high-contribution genes for further NMF.
#'
#' @param count A gene count matrix.
#' @param sgrna A cell-perturbation data frame. The colname of sgrna data frame must be 'cell' and 'perturbation'.
#' @param factor The number of factor in perturbmeta Default: 20.
#' @param nfeatures A parameter for extracting high-variable genes. Default: 3000.
#' @param npcs A number of PCs in PCA. Default: 20 (Same as number of factor).
#' @param filter The ratio of filtering cells that are perturbed non-significant. Default: 0.5.
#' @param ctrl The name of control perturbation. Default: NT.
#' @param VIP The threshold of extracting features from PLS-DA. Default: 1.
#' @param cluster The number of clusters that used for perturbmeta.
#' @param ncores The number of cores used for calculating NMF. Default: 0.5.
#'
#' @return A list contain the result of perturbmeta for each cluster
#' @export
#' @import Seurat
#' @import NMF
#' @import dplyr
#' @import doParallel
#'
#' @examples
#' res <- perturbmeta(count, sgrna,factor=20, nfeatures=3000, ctrl='NT', cluster=5)
perturbmeta <- function(count, sgrna, factor=20, nfeatures=3000, npcs=20, filter=0.5, ctrl='NT',VIP=1, cluster=5, ncores=0.5){

  if(sum(colnames(sgrna) %in% c('cell','perturbation')) != 2){
    stop('The sgrna data do not contain cell, perturbation. Please check sgrna data again!')
  }

  nam <- intersect(colnames(count), sgrna$cell)
  if(length(nam) == 0){
    stop('The count data and sgrna data contain different cells. Please check input data again!')
  }
  print(paste(paste('The count data and sgrna data have', length(nam), sep = ' '), 'cells with the same names',sep = ' '))

  sgrna <- sgrna[which(sgrna$cell %in% nam),]
  count <- count[,which(colnames(count) %in% nam)]

  rna <- CreateSeuratObject(counts = count)
  rna <- NormalizeData(rna) %>% FindVariableFeatures(nfeatures = features) %>% ScaleData(do.center = F)

  rna_filter <- filter_cell(counts=rna, sgrna=sgrna, filter=filter, npcs=npcs, nfeatures=nfeatures, ctrl=ctrl)
  print(paste(dim(rna_filter)[2], 'cells are retained for subsequent analysis',sep = ' '))

  perturbation_filter <- sgrna[which(sgrna$cell %in% colnames(rna_filter)),]
  pert<-as.character(unique(perturbation_filter$perturbation))

  nmf_data <- LayerData(rna_filter, assay = "RNA", layer = 'scale.data')

  num_cores <- detectCores()
  num_cores <- num_cores*0.5
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  rank_num = factor
  res_filter <- nmf(nmf_data, rank = rank_num, method="snmf/r", seed='nndsvd',.opt = 'vP')

  rna_filter@reductions$nmf@cell.embeddings<-t(coef(res_filter))
  rna_filter@reductions$nmf@feature.loadings<-basis(res_filter)
  rna_filter <- RunUMAP(rna_filter, reduction = 'nmf',dims = 1:rank_num,reduction.key = 'nmfumap',reduction.name = 'nmfumap') %>% FindNeighbors(reduction = 'nmf', dims = 1:rank_num) %>% FindClusters()

  idnt <- levels(rna_filter@active.ident)
  idnt <- idnt[1:cluster]
  variable_feature<-VariableFeatures(rna_filter)

  rna1<-CreateSeuratObject(counts = rna_filter@assays$RNA@counts)
  rna1@active.ident<-rna_filter@active.ident

  embedding_list <- list()
  loading_list <- list()
  feature_list <- list()

  for(i in 1:cluster){
    a<-subset(rna1, idents=idnt[i])
    count1<-data.frame(LayerData(a, assay = "RNA", layer = 'count'))
    count1<-count1[which(rownames(count1) %in% variable_feature),]

    plsda_res <- RunPLSDA(count=count1, sgrna=perturbation_filter, VIP=VIP, ncomp=npcs)
    embedding_list[[i]] <- plsda_res[[1]]
    loading_list[[i]] <- plsda_res[[2]]
    feature_list[[i]] <- plsda_res[[3]]
  }

  res_nmf <- list()
  for(i in 1:cluster){
    a<-subset(rna1, idents=idnt[i])#做了NMF
    nmf_data1<-nmf_data[,which(colnames(nmf_data) %in% colnames(a))]#不做PLS-DA
    nmf_data2<-nmf_data1[which(rownames(nmf_data1) %in% c(pert, feature_list[[i]][,1])),]#做了PLS-DA
    nmf_data2<-nmf_data2[which(rowSums(nmf_data2) > 0),]
    rank_num = factor#设置factor数量为20
    res1 <- nmf(nmf_data2, rank = rank_num, method="snmf/r", seed='nndsvd',.opt = 'vP')
    res_nmf[[i]]<-res1
    gc()
  }
  return(res_nmf)
}
