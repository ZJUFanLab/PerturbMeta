#' Plot the result of perturbmeta score
#' @description
#' Use dotplot to plot the result of perturbmeta score. It can plot the significant perturbmeta score (p < 0.05) for all cluster
#' by setting the plot parameter  to all. While it can plot all perturbmeta score from interested cluster by setting the plot
#' parameter to the number of interested cluster.
#'
#' @param score The first result of function perturbmeta_score.
#' @param p The second result of function perturbmeta_score.
#' @param cluster The number of clusters that used for perturbmeta.
#' @param plot A parameter for which cluster to plot dotplot. Default: 'all'.
#' @param nfact The number of factor in NMF. Default: 20.
#'
#' @return A dotplot for oerturbmeta score from interested cluster.
#' @export
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import tibble
#'
#' @examples
#' p<-perturbmeta_score_plot(res_score[[1]], res_score[[2]])
perturbmeta_score_plot <- function(score, p, cluster=5, plot='all', nfact=20){
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

  if(plot=='all'){
    data_all <- list()
    for (i in 1:cluster){
      data_all[[i]] <- PerturbMetascore_Bonferroni_score_list[[i]] %>% mutate(Group = paste0("Cluster",i-1,"_",factor))
    }
    data_all <- do.call("rbind",data_all)
    data_all_filter <- data_all[data_all$Bonferroni < 0.05,]

    ## 绘图
    data <- data_all_filter
    data$Bonferroni <- round(data$Bonferroni,3)
    Cluster_num = cluster
    p = ggplot(data,aes(Group,perturbation,size = -log10(Bonferroni) ))+
      geom_point(shape=21,aes(fill= PerturbMeta.score),position =position_dodge(0))+ ##shape = 20是没有包边的圆，21是有黑色包边的圆
      theme_minimal()+xlab(NULL)+ylab(NULL) +theme_bw()+
      scale_size_continuous(range=c(4,6))+theme_bw()+
      scale_fill_gradient2(low = "blue",mid ="white", midpoint = 0, high = "red")+
      theme(legend.position = "right",legend.box = "vertical", #图例位置
            legend.margin=margin(t= 0, unit='cm'),
            legend.spacing = unit(0,"in"),
            # 移除背景格子线
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            ##
            axis.text.x  = element_text(color="black",size=10,angle = 45,
                                        vjust = 1,hjust = 1),#x轴
            axis.text.y  = element_text(color="black",size=10),#y轴
            legend.text = element_text(size =10,color="black"),#图例
            legend.title = element_text(size =10,color="black", face = "bold"),#图例
            axis.title.y=element_text(vjust=1,size=10),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5) # 设置标题加粗并居中
      )+labs(x=" ",y = "",size = "Bonferroni adjusted \n p-value (-log10)")+
      ggtitle(paste0(Cluster_num, " clusters"))
    p
    return(p)
  }
  else{
    i = as.numeric(plot)
    data = PerturbMetascore_Bonferroni_score_list[[i]]
    data$factor <- gsub("factor","Factor",data$factor)
    data$factor <- factor(data$factor, levels = paste0("Factor",rep(1:nfact)))
    data$Bonferroni <- round(data$Bonferroni,3)
    data <- data %>% mutate(size_var = factor(ifelse(Bonferroni > 0.25, 1, ifelse(Bonferroni > 0.05, 2,3))))
    p = ggplot(data,aes(factor,perturbation,size = size_var ))+
      geom_point(shape=21,aes(fill= PerturbMeta.score),position =position_dodge(0))+ ##shape = 20是没有包边的圆，21是有黑色包边的圆
      theme_minimal()+xlab(NULL)+ylab(NULL) +theme_bw()+
      # scale_size_continuous(range=c(1,7.5))+theme_bw()+
      scale_fill_gradient2(low = "blue",mid ="white", midpoint = 0, high = "red")+
      theme(legend.position = "right",legend.box = "vertical", #图例位置
            legend.margin=margin(t= 0, unit='cm'),
            legend.spacing = unit(0,"in"),
            # 移除背景格子线
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            ##
            axis.text.x  = element_text(color="black",size=10,angle = 45,
                                        vjust = 1,hjust = 1),#x轴
            axis.text.y  = element_text(color="black",size=10),#y轴
            legend.text = element_text(size =10,color="black"),#图例
            legend.title = element_text(size =10,color="black", face = "bold"),#图例
            axis.title.y=element_text(vjust=1,size=10),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5) # 设置标题加粗并居中
      )+labs(x=" ",y = "",size = "Bonferroni adjusted \n p-value")+
      ggtitle(paste0("cluster",i-1,"")) +
      # 手动设置size的图例
      scale_size_manual(values = c(1, 3, 5), # 假设1, 2, 3是size_var的三个级别对应的值
                        breaks = c(1, 2, 3),
                        labels = c(">0.25", "0.05-0.25", "0-0.05")) # 设置圆为实心黑色
    p
    return(p)
  }
}
