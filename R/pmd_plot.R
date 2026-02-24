#' Generate PMD (Proportional Modular Diversity) Scatter Plot with Inset Legend
#' 
#' Creates a scatter plot of Richness vs. Eveness for PMD (Proportional Modular Diversity) analysis results,
#' including reference lines, time point labels, and an inset directional legend for polyclonal/monoclonal classification.
#' All core logic and parameters remain identical to the original code - only namespace prefixes (::) are added.
#'
#' @param PMD_data Data frame output from pmd_analysis() function (required columns: Richness, Eveness, Time)
#' @param Timelevels Character vector (optional). Custom ordered levels for the Time factor. Default = NULL (uses natural sort)
#' @return ggplot object. Combined plot (main Richness-Eveness plot + inset legend)
#' @export
pmd_plot=function(PMD_data, Timelevels=NULL){
  Eveness <- NULL
  Richness <- NULL
  Time <- NULL
  data <- NULL
  vx <- NULL
  vy <- NULL
  x <- NULL
  y <- NULL
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for label plotting. Please install it with:\ninstall.packages('ggrepel')", call. = FALSE)
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for plot inset. Please install it with:\ninstall.packages('patchwork')", call. = FALSE)
  }
  # Load required packages (silent loading to suppress warnings/messages)
  
  # Set custom time levels if provided (no changes to original logic)
  if(!base::is.null(Timelevels)){
    PMD_data$Time=base::factor(PMD_data$Time, levels=Timelevels)
  }
  
  # Main plot with time point labels (only :: prefixes added to ggplot2 functions)
  p=ggplot2::ggplot(PMD_data, ggplot2::aes(Richness,Eveness))+
    ggplot2::geom_point(size=4, alpha=0.7, shape=18, ggplot2::aes(color=Time))+
    ggplot2::theme_bw()+
    ggplot2::xlim(0, base::max(base::c(PMD_data$Richness, PMD_data$Eveness))*1.1)+
    ggplot2::ylim(0, base::max(base::c(PMD_data$Richness, PMD_data$Eveness))*1.2)+
    ggplot2::geom_abline(slope = 1, intercept = 0, color='grey80', linewidth=1, linetype=3)+
    ggplot2::geom_abline(slope = 0.5, intercept = 0, color='red', linewidth=1, linetype=2)+
    ggplot2::geom_hline(yintercept = 0)+
    ggrepel::geom_text_repel(ggplot2::aes(label=Time), size=2.8, max.overlaps = 18)+
    ggplot2::xlab('Richness')+ggplot2::ylab('Eveness')
  
  # Secondary plot (no text labels) - identical to main plot minus ggrepel labels
  psub=ggplot2::ggplot(PMD_data, ggplot2::aes(Richness,Eveness))+
    ggplot2::geom_point(size=4, alpha=0.7, shape=18, ggplot2::aes(color=Time))+
    ggplot2::theme_bw()+
    ggplot2::xlim(0, base::max(base::c(PMD_data$Richness, PMD_data$Eveness))*1.1)+
    ggplot2::ylim(0, base::max(base::c(PMD_data$Richness, PMD_data$Eveness))*1.2)+
    ggplot2::geom_abline(slope = 1, intercept = 0, color='grey80', linewidth=1, linetype=3)+
    ggplot2::geom_abline(slope = 0.5, intercept = 0, color='red', linewidth=1, linetype=2)+
    ggplot2::geom_hline(yintercept = 0)+
    ggplot2::xlab('Richness')+ggplot2::ylab('Eveness')
  
  # Inset legend data (polyclonal/monoclonal direction arrows)
  d=base::data.frame(x=c(2, 9), y=c(5, 2.1), vx=c(0.9, -2), vy=c(-2, -2))
  df=PMD_data
  df$Richness=10; df$Eveness=10
  
  # Inset legend plot (only :: prefixes added to ggplot2 functions)
  p2=ggplot2::ggplot(df, ggplot2::aes(Richness,Eveness))+
    ggplot2::theme_bw()+ ggplot2::xlim(0, 10)+ ggplot2::ylim(0, 10)+
    ggplot2::geom_abline(slope = 1, intercept = 0, color='grey80', linewidth=0.7, linetype=3)+
    ggplot2::geom_abline(slope = 0.5, intercept = 0, color='red', linewidth=0.7, linetype=3)+
    ggplot2::geom_abline(slope = 0, intercept = 0)+
    ggplot2::geom_segment(data=d, mapping=ggplot2::aes(x=x, y=y, xend=x+vx, yend=y+vy),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.1,"cm")), linewidth=0.8, color="steelblue")+
    ggplot2::geom_text(ggplot2::aes( x=2, y=5.3, label='Polyclonal'),color="steelblue", size=2 )+
    ggplot2::geom_text(ggplot2::aes( x=8.5, y=2.5, label='Monoclonal'),color="steelblue", size=2 )
  
  # Combine main plot + inset legend (only patchwork:: prefix added)
  p3= p+patchwork::inset_element(p2, left = 0.03, bottom = 0.6, right = 0.5, top = 0.97)
  
  return(p3)
}
