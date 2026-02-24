#' Plot Cancer-associated gene clone contribution
#' 
#' This function filters integration site data for cancer-associated genes (within specified distance/threshold)
#' and generates a dot plot of clone contribution percentages for these genes.
#' 
#' @param IS_raw Data frame with annotated integration site data (columns: nearest_gene_name, nearest_distance, Clone_contribution, Sample required)
#' @param Distance Numeric, maximum distance to cancer gene (default: 100000 bp)
#' @param threashold Numeric, minimum clone contribution threshold (default: 0.001)
#' @return ggplot object (dot plot of clone contribution for cancer-associated genes)
#' @export
is_in_CG_gene=function(IS_raw, Distance=100000, threashold=0.001){
  nearest_gene_name <- NULL
  Clone_contribution <- NULL
  Sample <- NULL
  index <- NULL
  cum_sum <- NULL
  Group <- NULL
  Freq <- NULL
  variable <- NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to generate AE gene plots. Install with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required to generate AE gene plots. Install with:\ninstall.packages('RColorBrewer')", call. = FALSE)
  }
  
  CG_gene=EPS$cancer_gene
  
  CG_sheet=IS_raw[IS_raw$nearest_gene_name %in% CG_gene, ]
  
  ############################################
  filter_gene_by_pos=function(gene_sheet){
    gene_sheet=gene_sheet[base::abs(gene_sheet$nearest_distance)<Distance,]
    return(gene_sheet)
  }
  
  CG_sheet=filter_gene_by_pos(CG_sheet)
  CG_sheet=CG_sheet[CG_sheet$Clone_contribution>=threashold,]
  
  CG_gene_sum=base::nrow(CG_sheet)
  if (CG_gene_sum==0){
    p=ggplot2::ggplot()+
      ggplot2::geom_abline(intercept = 0)+
      ggplot2::annotate("text", x=0, y=0, label= "No cancer gene")+  
      ggplot2::theme_classic()
    return(p)
  }
  
  natural_sort <- function(x) {
    numeric_part <- base::as.numeric(base::gsub("[^0-9.]", "", x))
    unit_part <- base::gsub("[0-9.]", "", x)
    base::order(numeric_part, unit_part)
  }
  
  CG_sheet$Sample=base::factor(CG_sheet$Sample, 
                               levels=base::unique(IS_raw$Sample)[natural_sort(base::unique(IS_raw$Sample))])
  
  ################################################
  mycol=base::c(RColorBrewer::brewer.pal(n = 9,name = 'Set1'),
                RColorBrewer::brewer.pal(n = 8,name = 'Paired'),
                RColorBrewer::brewer.pal(n = 8,name = 'Dark2'))
  
  if (base::length(base::unique(CG_sheet$Sample))>20){
    p=  ggplot2::ggplot(CG_sheet, ggplot2::aes(nearest_gene_name, Clone_contribution*100))+
      ggplot2::geom_point(ggplot2::aes(nearest_gene_name, Clone_contribution*100), size=5, alpha=0.6, shape=18)+
      ggplot2::theme_bw()+
      ggplot2::ylim(0, base::max(CG_sheet$Clone_contribution*100*1.1, 5, na.rm = TRUE))+
      ggplot2::scale_color_manual(values = mycol)+
      ggplot2::xlab('Cancer associated genes')+ggplot2::ylab('Clone contribution [%]')+
      ggplot2::theme(legend.position = "right")+
      ggplot2::labs(color='')+
      ggplot2::coord_flip()
  }
  
  if (base::length(base::unique(CG_sheet$Sample))>3 & base::length(base::unique(CG_sheet$Sample))<=20){
    p= ggplot2::ggplot(CG_sheet, ggplot2::aes(nearest_gene_name, Clone_contribution*100))+
      ggplot2::geom_point(ggplot2::aes(nearest_gene_name, Clone_contribution*100, color=Sample), size=5, alpha=0.6, shape=18)+
      ggplot2::theme_bw()+
      ggplot2::ylim(0, base::max(CG_sheet$Clone_contribution*100*1.1, 5, na.rm = TRUE))+
      ggplot2::scale_color_manual(values = mycol)+
      ggplot2::xlab('Cancer associated genes')+ggplot2::ylab('Clone contribution [%]')+
      ggplot2::theme(legend.position = "right")+
      ggplot2::labs(color='')+
      ggplot2::coord_flip()
  } else if(base::length(base::unique(CG_sheet$Sample))<=3){
    p= ggplot2::ggplot(CG_sheet, ggplot2::aes(nearest_gene_name, Clone_contribution*100))+
      ggplot2::geom_point(ggplot2::aes(nearest_gene_name, Clone_contribution*100, color=Sample), size=5, alpha=0.6, shape=18)+
      ggplot2::theme_bw()+
      ggplot2::ylim(0, base::max(CG_sheet$Clone_contribution*100*1.1, 5, na.rm = TRUE))+
      ggplot2::scale_color_manual(values = mycol)+
      ggplot2::xlab('Cancer associated genes')+ggplot2::ylab('Clone contribution [%]')+
      ggplot2::theme(legend.position = "top")+
      ggplot2::labs(color='')+
      ggplot2::coord_flip()
  }
  
  return(p)
}


