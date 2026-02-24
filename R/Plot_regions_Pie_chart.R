#' Plot Region-wise Donut Charts
#' @param Region_data Named list of data frames with Product/Share/Percentage/Time columns
#' @param Timelevels Character vector to subset time levels (optional)
#' @return Arranged ggplot object of donut charts
#' @export
plot_regions <- function(Region_data, Timelevels = NULL) {
  
  ggdonutchart <- NULL
  ggarrange <- NULL
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required. Install with: install.packages('ggpubr')", call. = FALSE)
  }
  
  natural_sort <- function(x) {
    numeric_part <- as.numeric(gsub("[^0-9.]", "", x))
    unit_part <- gsub("[0-9.]", "", x)
    order(numeric_part, unit_part)
  }
  
  Region_data = Region_data[natural_sort(names(Region_data))]
  if(!is.null(Timelevels)){
    Region_data = Region_data[Timelevels]
  }
  
  # Custom color scheme
  colors <- c(
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FECA57",
    "#FF9FF3", "#54A0FF", "#5F27CD", "#66C2A5", "#FC8D62"
  )
  
  for(i in names(Region_data)){
    Region_data[[i]]$angle <- 90 - 360 * cumsum(c(0, Region_data[[i]]$Share[-nrow(Region_data[[i]])])) / sum(Region_data[[i]]$Share)
    Region_data[[i]]$hjust <- ifelse(Region_data[[i]]$angle < -90, 1, 0)
    Region_data[[i]]$angle <- ifelse(Region_data[[i]]$angle < -90, Region_data[[i]]$angle + 180, Region_data[[i]]$angle)
    Region_data[[i]]$Order=1:nrow(Region_data[[i]])
    Region_data[[i]]$name=paste0(Region_data[[i]]$Product, '\n', Region_data[[i]]$Percentage,'%')
  }
  
  # Donut chart(ggpubr::ggdonutchart)
  my_ggdounut=list()
  for(i in names(Region_data)){
    my_ggdounut[[i]]=ggpubr::ggdonutchart(
      data = Region_data[[i]][Region_data[[i]]$Percentage>0,],
      x = "Share", fill = "Product",
      label = rep('', nrow(Region_data[[i]][Region_data[[i]]$Percentage>0,])),
      color = "white",  palette = 'Set1',
      lab.pos = c('out'),  lab.adjust = 0.04, lab.font = 3,
      title = Region_data[[i]]$Time[1]
    )
  }
  
  # Arrange plots（ggpubr::ggarrange）
  if (length(my_ggdounut)==2){
    p=ggpubr::ggarrange(my_ggdounut[[1]],my_ggdounut[[2]] ,
                        ncol = 2, nrow = 1, align = "h",  common.legend =TRUE)
  }
  if(length(my_ggdounut)==3){
    p=ggpubr::ggarrange(my_ggdounut[[1]],my_ggdounut[[2]] ,my_ggdounut[[3]],
                        ncol = 3, nrow = 1, align = "h",  common.legend =TRUE)
  }
  if(length(my_ggdounut)==4){
    p=ggpubr::ggarrange(my_ggdounut[[1]],my_ggdounut[[2]] ,my_ggdounut[[3]],
                        my_ggdounut[[4]],
                        ncol = 2, nrow = 2, align = "h",  common.legend =TRUE)
  }
  if(length(my_ggdounut)==5){
    p=ggpubr::ggarrange(my_ggdounut[[1]],my_ggdounut[[2]] ,my_ggdounut[[3]],
                        my_ggdounut[[4]],my_ggdounut[[5]],
                        ncol = 2, nrow = 2, align = "h",  common.legend =TRUE)
  }
  if(length(my_ggdounut)==6){
    p=ggpubr::ggarrange(my_ggdounut[[1]],my_ggdounut[[2]] ,my_ggdounut[[3]],
                        my_ggdounut[[4]],my_ggdounut[[5]],my_ggdounut[[6]],
                        ncol = 3, nrow = 2, align = "h",  common.legend =TRUE)
  }
  return(p)
}
