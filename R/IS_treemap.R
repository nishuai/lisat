#' Generate treemap of integration site clone contribution
#' 
#' This function creates a treemap visualization of the top 1000 integration site (IS) clone contributions,
#' grouped by patient time points with custom color perturbation.
#' 
#' @param IS_raw Data frame containing IS data (columns: Sample, Locus, Clone_contribution required)
#' @param Patient_timepoint Data frame mapping Sample_ID to Time_Point (columns: Sample_ID, Time_Point required)
#' @param Timelevels Character vector, optional custom order of time points (default: NULL, natural sort)
#' @return ggplot object (treemap of IS clone contributions)
#' @export
IS_treemap=function(IS_raw=IS_raw, Patient_timepoint = Patient_timepoint, Timelevels=NULL){
  Clone_contribution <- NULL
  Locus <- NULL
  Sample <- NULL
  # ====================== 第一步：加载依赖包（确保函数可用） ======================
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package 'colorspace' is required. Please reinstall the package.", call. = FALSE)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please reinstall the package.", call. = FALSE)
  }
  
  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop("Package 'treemapify' is required. Please reinstall the package.", call. = FALSE)
  }
  
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("Package 'viridisLite' is required. Please reinstall the package.", call. = FALSE)
  }
  
  # base_colors= c("#E0F7FA", "#81D4FA", "#4A69A1", "#01579B","#4FC3F7","#1E8AFE")
  # base_colors=base::rep(base_colors, 2000)
  # base::set.seed(100)
  # 
  # res <- base::character()
  # rgb_vals <- grDevices::col2rgb(base_colors) / 255  # 得到的是 3x1 矩阵
  # sRGB_obj <- colorspace::sRGB(rgb_vals[1, ], rgb_vals[2, ], rgb_vals[3, ])
  # 
  # hcl_col <- methods::as(sRGB_obj, "polarLUV")
  # 
  # H <- hcl_col @ coords[, "H"]
  # C <- hcl_col @ coords[, "C"]
  # L <- hcl_col @ coords[, "L"]
  # 
  # new_H <- H + stats::runif(base::length(H), -5, 5)  # 色相扰动 ±15
  # new_C <- C + stats::runif(base::length(C), -5, 5)  # 饱和度扰动 ±15
  # new_L <- L + stats::runif(base::length(L), -5, 5)  # 亮度扰动 ±15
  # 
  # new_col <- colorspace::polarLUV(L = new_L, C = new_C, H = new_H)
  # res <- base::c(res, colorspace::hex(new_col))
  # 
  # res=res[!base::is.na(res)]
  res=EPS$res
  
  # ------------------------------------------
  base::set.seed(28)
  IS_raw=IS_raw[,base::c('Sample','Locus','Clone_contribution')]
  IS_raw$Locus=base::paste0(IS_raw$Locus, base::seq_len(base::nrow(IS_raw)))
  
  groups <- base::unique(IS_raw$Sample)
  for(i in groups){
    IS_raw$Sample[IS_raw$Sample==i] =
      Patient_timepoint$Time_Point[Patient_timepoint$Sample_ID==i]
  }
  
  natural_sort <- function(x) {
    numeric_part <- base::as.numeric(base::gsub("[^0-9.]", "", x))
    unit_part <- base::gsub("[0-9.]", "", x)
    base::order(numeric_part, unit_part)
  }
  
  groups <- base::unique(IS_raw$Sample)
  if(base::is.null(Timelevels)){
    groups=groups[natural_sort(groups)]
  } else{
    groups=Timelevels
  }
  IS_raw$Sample=base::factor(IS_raw$Sample, levels=groups)
  
  IS_raw=base::split(IS_raw, f=IS_raw$Sample)
  IS_raw=base::lapply(IS_raw, function(x) x[1:1000,])
  IS_raw=base::do.call(base::rbind, IS_raw)
  IS_raw=IS_raw[!base::is.na(IS_raw$Sample),]
  
  p=ggplot2::ggplot(IS_raw, ggplot2::aes(
    area = Clone_contribution, 
    fill = Locus,          
    subgroup = Sample
  )) +
    treemapify::geom_treemap(radius =grid::unit(0.2, "npc"),
                             color = "white",    
                             size = 0.5
    )+
    ggplot2::scale_fill_manual(values = res) +
    ggplot2::facet_wrap(~ Sample, ncol = base::length(base::unique(IS_raw$Sample))) +
    ggplot2::labs(
      title = "Top 1000 UIS Relative Percentages",
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::coord_fixed()  
  
  return(p)
}
