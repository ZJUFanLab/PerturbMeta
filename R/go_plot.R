#' Plot the result of GO enrichment analysis for interested cluster
#' @description
#' Use dotplot to plot the result of GO enrichment analysis. It can plot the GO enrichment results (p < 0.05) for interested cluster
#' by setting the plot parameter  to the number of interested cluster.
#'
#' @param go_res The result from function enrich_go.
#' @param score The first result of function perturbmeta_score.
#' @param p The second result of function perturbmeta_score .
#' @param cluster The number of clusters that used for perturbmeta.
#' @param top The number of top GO pathways. Default: 1.
#' @param factor The number of factor in perturbmeta Default: 20.
#' @param plot A parameter for which cluster to plot dotplot.
#'
#' @return A dotplot for GO enrichment analysis for interested cluster
#' @export
#' @import tidyr
#' @import tibble
#' @import ggpubr
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' p<-go_plot(go_res,res_score[[1]],res_score[[2]],top=5,plot=1)
go_plot <- function(go_res, score, p, cluster=5, top=5, factor=20, plot=1){

  Allclusters_GO_df <-list()
  for (i in 1:cluster){
    temp_GO <- go_res[[i]]
    for (j in 1:factor){
      temp_GO[[j]] <- temp_GO[[j]] %>% mutate(Factor = paste0("factor",j))
    }
    temp_GO <- do.call("rbind",temp_GO)
    Allclusters_GO_df[[i]] <- temp_GO %>% mutate(Cluster = paste0(i))
  }
  Allclusters_GO_df <- do.call("rbind", Allclusters_GO_df) ##合并一个完整数据框

  Allclusters_GO_df$Group <- with(Allclusters_GO_df,
                                  paste0(Cluster,"_", Factor))

  PerturbMetascore_Bonferroni_score_list <- list()
  for (i in seq_along(score)){
    PerturbMetascore1 <- score[[i]] %>%
      data.frame() %>%
      rownames_to_column(var = "perturbation") %>%
      pivot_longer(cols = contains("factor"),names_to = "factor",values_to = "PerturbMeta score")
    unique_idx1 <- paste0(PerturbMetascore1$perturbation, PerturbMetascore1$factor)

    Bonferroni_score1 <- p[[i]] %>%
      t() %>%
      data.frame() %>%
      rownames_to_column(var = "perturbation") %>%
      pivot_longer(cols = contains("factor"),names_to = "factor",values_to = "Bonferroni")
    unique_idx2 <- paste0(Bonferroni_score1$perturbation, Bonferroni_score1$factor)

    PerturbMetascore_Bonferroni_score <- PerturbMetascore1 %>% mutate("Bonferroni" = Bonferroni_score1$Bonferroni[match(unique_idx1, unique_idx2)]) %>% data.frame()
    PerturbMetascore_Bonferroni_score_list[[i]] <- PerturbMetascore_Bonferroni_score
  }

  clusters_GO <- list()
  for (i in 1:5){
    tmp_df <- PerturbMetascore_Bonferroni_score_list[[i]]
    clu <- paste0("Cluster",i-1)
    tmp_df_sig <- tmp_df %>% filter(Bonferroni < 0.05)
    Group_id <- paste0(cluster,"_",tmp_df_sig$factor)
    GO_df <- Allclusters_GO_df[Allclusters_GO_df$Group %in% Group_id,]
    clusters_GO[[clu]] <- GO_df %>% filter(pvalue < 0.05)
  }

  topn = top ##绘制top GO terms的数量
  i=plot ##指定哪个Cluster, i = 1 表示指定了cluster0
  Cluster = names(clusters_GO)[i]
  Factor_GO_df <- clusters_GO[[Cluster]]
  Factor_GO_df$Factor <- gsub("factor","Factor_",Factor_GO_df$Factor)
  Factor_GO_df_topn <- Factor_GO_df %>% group_by(Factor) %>% arrange(Factor, pvalue) %>% slice_head(n = topn)
  Factor_GO_df_topn$Description <- factor(Factor_GO_df_topn$Description,levels = unique(Factor_GO_df_topn$Description))


  ## 绘图
  p = ggplot(Factor_GO_df_topn,aes(Factor,Description,size = Count))+
    geom_point(shape=21,aes(fill= (-log10(pvalue))),position =position_dodge(0))+ ##shape = 20是没有包边的圆，21是有黑色包边的圆
    theme_minimal()+xlab(NULL)+ylab(NULL) +
    scale_size_continuous(range=c(4,8))+
    theme_bw()+
    scale_fill_gradient2(low = "#FDE300", high = "#002060")+
    theme(legend.position = "right",legend.box = "vertical", #图例位置
          legend.margin=margin(t= 0, unit='cm'),
          legend.spacing = unit(0,"in"),
          # 移除背景格子线
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          ##
          axis.text.x  = element_text(color="black",size=14,angle = 45,
                                      vjust = 1,hjust = 1),#x轴
          axis.text.y  = element_text(color="black",size=14),#y轴
          legend.text = element_text(size =14,color="black"),#图例
          legend.title = element_text(size =14,color="black"),#图例
          axis.title.y=element_text(vjust=1,size=14),
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5) # 设置标题加粗并居中
    )+labs(x=" ",y = "",fill = expression(paste("-", "log"[10], "(P value)")), size = "Enrichment ratio")+
    ggtitle(paste0(Cluster,'',sep=''))
  p
  return(p)
}
