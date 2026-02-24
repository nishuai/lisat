#' Plot Richness & Evenness Dual Y-Axis Line Chart
#' 
#' Creates a polished dual Y-axis line chart to visualize clonal richness and evenness over time,
#' with automatic scaling between axes, customizable styling, and optional data labels.
#' All core functionality and parameters remain identical to the original code - only namespace prefixes (::) added.
#'
#' @param PMD_data Data frame containing time, richness, and evenness data (required columns specified by time_col/richness_col/evenness_col)
#' @param time_col Character (default = "Time"). Name of column containing time points.
#' @param richness_col Character (default = "Richness"). Name of column containing richness values.
#' @param evenness_col Character (default = "Eveness"). Name of column containing evenness values (note: intentional spelling match to original code).
#' @param plot_title Character (default = "Clonal eveness over time"). Main plot title (spelling preserved as original).
#' @param subtitle Character (optional). Plot subtitle (default = NULL).
#' @param richness_color Character (default = "#3366CC"). Hex color code for richness line/points/labels.
#' @param evenness_color Character (default = "#CC6677"). Hex color code for evenness line/points/labels.
#' @param show_labels Logical (default = TRUE). Whether to display numeric labels on data points.
#' @param Timelevels Character vector (optional). Custom ordered levels for time factor (overrides default ordering).
#' @return ggplot object. Dual Y-axis line chart of richness (primary) and evenness (secondary) over time.
#' @export
plot_richness_evenness <- function(PMD_data, time_col = "Time", richness_col = "Richness",
                                   evenness_col = "Eveness", plot_title = "Clonal eveness over time",
                                   subtitle = NULL, richness_color = "#3366CC", evenness_color = "#CC6677",
                                   show_labels = TRUE, Timelevels=NULL) {
  .data=NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for scale adjustment. Please install it with:\ninstall.packages('scales')", call. = FALSE)
  }
  
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for color palette. Please install it with:\ninstall.packages('RColorBrewer')", call. = FALSE)
  }
  # Load required packages (preserved original package loading logic)
  
  # Check for required columns in input dataframe
  if (!base::all(base::c(time_col, richness_col, evenness_col) %in% base::colnames(PMD_data))) {
    stop(base::paste("Dataframe must contain", time_col, richness_col, "and", evenness_col, "column"))
  }
  
  # Ensure time column is ordered factor (preserved original logic)
  if (!base::is.factor(PMD_data[[time_col]])) {
    PMD_data[[time_col]] <- base::factor(PMD_data[[time_col]], levels = base::unique(PMD_data[[time_col]]))
  }
  
  if(!base::is.null(Timelevels)){
    PMD_data[[time_col]] <- base::factor(PMD_data[[time_col]], levels = Timelevels)
  }
  
  # Set modern minimal theme (only :: prefixes added to ggplot2 functions)
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 14) +
                       ggplot2::theme(
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.grid.major.x = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(color = "gray30"),
                         axis.ticks = ggplot2::element_line(color = "gray30"),
                         legend.position = "top",
                         legend.title = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40"),
                         plot.caption = ggplot2::element_text(color = "gray50")
                       ))
  
  # Define line/point styling (preserved original values)
  line_size <- 1.2
  point_size <- 3.5
  point_shape <- 21  # Filled circle with border
  
  # Calculate scale factor to align y-axis ranges (preserved original logic)
  max_richness <- base::max(PMD_data[[richness_col]])
  max_evenness <- base::max(PMD_data[[evenness_col]])
  scale_factor <- max_richness / max_evenness
  
  # Create base plot (only :: prefixes added to ggplot2 functions)
   p=ggplot2::ggplot(PMD_data, ggplot2::aes(x = .data[[time_col]])) +
    
    # Add Richness line and points
    ggplot2::geom_line(ggplot2::aes(y = .data[[richness_col]], color = "'Richness'"), linewidth = line_size, group=1) +
    ggplot2::geom_point(ggplot2::aes(y = .data[[richness_col]] ),
                        shape = point_shape, size = point_size, stroke = 1.2) +
    
    # Add Evenness line (scaled) and points (preserved original scaling logic)
    ggplot2::geom_line(ggplot2::aes(y = .data[[evenness_col]]*scale_factor, 
                                    color = "'Evenness'"),
                       linewidth = line_size, linetype = "dashed",group=1) +
    ggplot2::geom_point(ggplot2::aes(y = .data[[evenness_col]]*scale_factor),
                        shape = point_shape, size = point_size, stroke = 1.2) +
    
    # Set y-axis limits (preserved original expansion logic)
    ggplot2::expand_limits(y = base::c(0, base::max(max_richness, max_evenness * scale_factor) * 1.1)) +
    
    # Add plot titles/labels (preserved original text)
    ggplot2::labs(
      title = base::paste0(plot_title ),
      subtitle = subtitle,
      x = "Time",
      caption = "Time post infusion"
    ) +
    
    # Define color/fill scales (preserved original values)
    ggplot2::scale_color_manual(
      values = base::c( "Evenness" = evenness_color, "Richness" = richness_color),
      labels = base::c("Evenness", "Richness")
    ) +
    ggplot2::scale_fill_manual(
      values = base::c("Richness" = "white", "Evenness" = "white"),
      labels = base::c("Richness", "Evenness")
    ) +
    
    # Add secondary Y-axis for Evenness (preserved original scaling)
    ggplot2::scale_y_continuous(
      name = "Richness",
      sec.axis = ggplot2::sec_axis(~ . / scale_factor, name = "Evenness")
    )
  
  # Add optional data labels (preserved original logic)
  if (show_labels) {
   p= p +
      # Richness value labels
      ggplot2::geom_text(ggplot2::aes(y = .data[[richness_col]],
                                      label = round(.data[[richness_col]], 2)),
                         vjust = -0.8, size = 4, fontface = "bold", color = richness_color) +
      # Evenness value labels (scaled back to original values)
      ggplot2::geom_text(ggplot2::aes(y = .data[[evenness_col]]*scale_factor,
                                      label = round(.data[[evenness_col]], 2)),
                         vjust = 2.5, size = 4, fontface = "bold", color = evenness_color)
  }
  
  return(p)
}
