#' GO enrichment analysis
#' @description
#' GO enrichment analysis was performed based on high-contributing genes from the perturbmeta results for interested cluster.
#'
#' @param nmf The result from function perturbmeta.
#' @param cluster The interested cluster for GO enrichment analysis. Default: 5.
#' @param factor The number of factor in NMF. Default: 20.
#' @param sp The species for gene count matrix. Defalut: org.Hs.eg.db.
#' @param top The number of high-contribution features extracted for GO enrichment analysis. Default: 50.
#'
#' @return A list contain the results of Go enrichment analysis for interested clusters.
#' @export
#' @importFrom clusterProfiler enrichGO
#' @importFrom AnnotationDbi mapIds
#' @importFrom DOSE setReadable
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom NMF extractFeatures
#'
#' @examples
#' go_res<-enrich_go(res_nmf, cluster=5, factor=20, sp='org.Hs.eg.db', top=50)
enrich_go <-function(nmf, cluster=5, factor=20, sp='org.Hs.eg.db', top=50){
  cluster_go_list <- list()
  for(i in 1:cluster){
    res<-nmf[[i]]
    f<-NMF::extractFeatures(res,top)
    f <- lapply(f, function(x) rownames(res)[x])
    f <- do.call('cbind',f)
    colnames(f)<-c(paste0('factor',1:factor))

    go_list<-list()

    if(sp=='org.Hs.eg.db'){
      for(j in 1:factor){
        Ensembl_ID <- f[,j]
        gene_entrezid<-AnnotationDbi::mapIds(x = org.Hs.eg.db,
                              keys = Ensembl_ID,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first")
        gene_entrezid<-gene_entrezid[!is.na(gene_entrezid)]
        go_enrich <- clusterProfiler::enrichGO(gene = gene_entrezid,
                              OrgDb = 'org.Hs.eg.db',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)
        go_enrich2<-DOSE::setReadable(go_enrich,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
        go_enrich2<-as.data.frame(go_enrich2)
        gc()
        go_list[[j]]<-go_enrich2
      }
    }
    else{
      for(j in 1:factor){
        Ensembl_ID <- f[,j]
        gene_entrezid<-AnnotationDbi::mapIds(x = org.Mm.eg.db,
                              keys = Ensembl_ID,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first")
        gene_entrezid<-gene_entrezid[!is.na(gene_entrezid)]
        go_enrich <- clusterProfiler::enrichGO(gene = gene_entrezid,
                              OrgDb = 'org.Mm.eg.db',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 0.5)
        go_enrich2<-DOSE::setReadable(go_enrich,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
        go_enrich2<-as.data.frame(go_enrich2)
        gc()
        go_list[[j]]<-go_enrich2
      }
    }
    cluster_go_list[[i]] <- go_list
    names(cluster_go_list)[i] <- paste('cluster',i-1,sep = '')
  }
  return(cluster_go_list)
}
