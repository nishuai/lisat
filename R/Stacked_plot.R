#' Generate Linked Timepoint Sankey + Stacked Bar Chart
#' 
#' Creates a highly customizable combined Sankey-flow + stacked bar chart to visualize clonal proportion changes across timepoints,
#' with manual control over flow polygon shapes and precise formatting of top integration sites (top 10 + "Others" category).
#' All core logic and data processing steps remain identical to the original code - only namespace prefixes (::) added and lag() fixed.
#'
#' @param IS_raw Data frame containing integration site data (required columns: Clone_contribution, Sample, nearest_gene_name, Chr, Locus)
#' @param Patient_timepoint Data frame mapping Sample_ID to Time_Point (columns: Sample_ID, Time_Point required)
#' @param Timelevels Character vector (optional). Custom ordered levels for time points (overrides natural sort). Default = NULL.
#' @return ggplot object. Combined Sankey-flow + stacked bar chart of top 10 integration site proportions across timepoints.
#' @importFrom dplyr group_by mutate ungroup filter arrange intersect summarise lag
#' @export
Linked_timepoints=function(IS_raw, Patient_timepoint, Timelevels=NULL){
  # Load required packages (silent loading to suppress warnings/messages)
  
  
  time_point <- value <- total <- percentage <- cum_percentage <- NULL
  x_pos <- category <- x <- y <- flow_id <- cum_percentage_lag <- NULL
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with: install.packages('dplyr')", call. = FALSE)
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required. Install with: install.packages('RColorBrewer')", call. = FALSE)
  }
  # Initialize list to store top 11 integration sites per sample
  top21_all=base::list()
  
  # Sort IS data by clone contribution (descending) and split by Sample
  IS_raw=IS_raw[base::order(IS_raw$Clone_contribution, decreasing = TRUE), ]
  IS_raw=base::split(IS_raw, f=IS_raw$Sample)
  sample_names=base::names(IS_raw)
  
  # Process each sample to extract top 10 + "Others" category
  for(i in sample_names){
    top21=IS_raw[[i]][1:11,]
    top21$Clone_contribution=top21$Clone_contribution*100
    
    # Calculate "Others" proportion (100% minus sum of top 10)
    top21$Clone_contribution[11]=100-base::sum(top21$Clone_contribution[1:10])
    
    # Handle NA values and standardize "Others" category
    top21$nearest_gene_name[base::is.na(top21$Clone_contribution)]='Others'
    top21$Clone_contribution[base::is.na(top21$Clone_contribution)]=0.001
    top21$ID=base::paste0(top21$nearest_gene_name,'_',top21$Chr,'_', top21$Locus)
    
    # Restructure data frame
    top21=top21[,base::c('Clone_contribution','ID')]
    top21$Sample_ID=i
    top21=top21[,base::c(2,3,1)]
    base::names(top21)=base::c('ID','Timepoint','Value')
    
    # Finalize "Others" category and remove duplicates
    top21$ID[11]='Others'
    top21$ID[base::grepl('Others', top21$ID)]='Others'
    top21=top21[!base::duplicated(top21$ID),]
    top21=top21[base::nrow(top21):1,]
    top21_all[[i]]=top21
  }
  
  # Combine all samples and map Sample_ID to Time_Point
  top21_all=base::do.call(base::rbind, top21_all)
  for(i in 1:base::nrow(top21_all)){
    if(base::sum(Patient_timepoint$Sample_ID==top21_all$Timepoint[i])==1){
      top21_all$Timepoint[i]=Patient_timepoint$Time_Point[Patient_timepoint$Sample_ID==top21_all$Timepoint[i]]
    }
  }
  
  # Add dummy timepoints for missing levels (to preserve Timelevels order in plot)
  if(base::length(base::unique(top21_all$Timepoint))< base::length(Timelevels)){
    missing_levels=Timelevels[! Timelevels %in% top21_all$Timepoint]
    for(i in missing_levels){
      top21_all=base::rbind(top21_all, top21_all[1,])
      top21_all$Timepoint[base::nrow(top21_all)]=i
      top21_all$Value[base::nrow(top21_all)]=0
    }
  }
  
  # Rename columns for consistency
  df_raw=top21_all
  base::names(df_raw)=base::c('category','time_point','value')
  
  # Order categories (custom order: top non-Others + Others last)
  category_names=base::unique(df_raw$category[base::order(df_raw$value, decreasing = TRUE)])
  category_names=category_names[category_names!='Others']
  category_names=base::c(category_names,'Others')
  df_raw$category=base::factor(df_raw$category, levels=category_names)
  
  #################### Set color palette (preserved original color logic) ####################
  colors <- base::c("grey60","#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4",
                    "#FECA57", "#FF9FF3", "#54A0FF", "#5F27CD")
  colors=base::c(colors, RColorBrewer::brewer.pal('Set2', n=8))
  colors=base::c(colors, RColorBrewer::brewer.pal('Set1', n=8))
  colors=base::c(colors, RColorBrewer::brewer.pal('RdBu', n=9))
  colors=colors[-(30:31)]
  colors=base::c(colors, RColorBrewer::brewer.pal('Dark2', n=8))
  colors=base::c(colors, RColorBrewer::brewer.pal('Accent', n=8))
  
  # Set random seed for reproducible color assignment
  base::set.seed(123)
  
  # Natural sort function for time points (preserved original logic)
  natural_sort <- function(x) {
    # Extract numeric part from time point labels
    numeric_part <- base::as.numeric(base::gsub("[^0-9.]", "", x))
    # Extract non-numeric (unit) part
    unit_part <- base::gsub("[0-9.]", "", x)
    # Order by numeric then unit part
    base::order(numeric_part, unit_part)
  }
  
  # Apply natural sort to time points (override with Timelevels if provided)
  df_raw$time_point=base::factor(df_raw$time_point,
                                 levels = base::unique(df_raw$time_point)[natural_sort(base::unique(df_raw$time_point))])
  if(!base::is.null(Timelevels)){
    df_raw$time_point=base::factor(df_raw$time_point, levels = Timelevels)
  }
  
  # Calculate percentages and cumulative values for stacked bars
  # CORE FIX: lag() is from dplyr, not tidyr - removed tidyr dependency
  df <- df_raw %>%
    dplyr::group_by(time_point) %>%
    dplyr::mutate(
      total = base::sum(value),
      percentage = value / total * 100,
      cum_percentage = base::cumsum(percentage),
      cum_percentage_lag = dplyr::lag(cum_percentage, default = 0)  # Fixed: dplyr::lag() instead of tidyr::lag()
    ) %>%
    dplyr::ungroup() %>%
    # Add numeric x positions for plotting (1-4 for time points)
    dplyr::mutate(x_pos = dplyr::case_when(
      time_point == base::levels(time_point)[1] ~ 1,
      time_point == base::levels(time_point)[2] ~ 2,
      time_point == base::levels(time_point)[3] ~ 3,
      time_point == base::levels(time_point)[4] ~ 4
    ))
  
  # Assign fixed colors to categories
  all_categories <- base::unique(df$category)
  base::names(colors) <- all_categories
  
  #################### Core Function: Create Sankey Flow Polygons ####################
  create_flow_polygons <- function(df) {
    x_pos <- category <- NULL
    flow_polygons <- base::data.frame()
    time_levels <- base::sort(base::unique(df$x_pos))
    
    # Iterate over adjacent time point pairs
    for(i in 1:(base::length(time_levels)-1)) {
      current_time <- time_levels[i]
      next_time <- time_levels[i+1]
      
      # Get data for current and next time points
      current_data <- df %>% dplyr::filter(x_pos == current_time)
      next_data <- df %>% dplyr::filter(x_pos == next_time)
      
      # Find common categories between time points
      common_categories <- dplyr::intersect(current_data$category, next_data$category)
      
      # Create flow polygons for each common category
      for(cat in common_categories) {
        current_cat <- current_data %>% dplyr::filter(category == cat)
        next_cat <- next_data %>% dplyr::filter(category == cat)
        
        if(base::nrow(current_cat) > 0 && base::nrow(next_cat) > 0) {
          # Get y-coordinates (bottom/top boundaries) for current time point
          y1_bottom <- current_cat$cum_percentage_lag
          y1_top <- current_cat$cum_percentage
          
          # Get y-coordinates for next time point
          y2_bottom <- next_cat$cum_percentage_lag
          y2_top <- next_cat$cum_percentage
          
          # Create 5 vertices for flow polygon (closed quadrilateral)
          polygon_x <- base::c(current_time + 0.15,     # Bottom-left (current time right edge)
                               next_time - 0.15,        # Bottom-right (next time left edge)
                               next_time - 0.15,        # Top-right
                               current_time + 0.15,     # Top-left
                               current_time + 0.15)     # Close polygon
          
          polygon_y <- base::c(y1_bottom,    # Bottom-left
                               y2_bottom,     # Bottom-right
                               y2_top,        # Top-right
                               y1_top,        # Top-left
                               y1_bottom)     # Close polygon
          
          # Create polygon data frame
          polygon_data <- base::data.frame(
            x = polygon_x,
            y = polygon_y,
            category = cat,
            flow_id = base::paste(cat, current_time, next_time, sep = "_"),
            stringsAsFactors = FALSE
          )
          
          flow_polygons <- base::rbind(flow_polygons, polygon_data)
        }
      }
    }
    
    return(flow_polygons)
  }
  
  # Generate flow polygon data for Sankey connections
  flow_data <- create_flow_polygons(df)
  
  #################### Create Main Plot (Sankey + Stacked Bars) ####################
  p <- ggplot2::ggplot() +
    # 1. Plot Sankey flow polygons (bottom layer)
    ggplot2::geom_polygon(data = flow_data,
                          ggplot2::aes(x = x, y = y, fill = category, group = flow_id),
                          alpha = 0.6,  # Transparency to show flow under bars
                          color = NA) +  # Remove polygon borders
    
    # 2. Plot stacked bar chart (top layer)
    ggplot2::geom_rect(data = df,
                       ggplot2::aes(xmin = x_pos - 0.15, xmax = x_pos + 0.15,
                                    ymin = cum_percentage_lag, ymax = cum_percentage,
                                    fill = category),
                       linewidth = 0.5,          # Border thickness
                       alpha = 0.9) +
    
    # 3. Set color palette (preserved original color assignment)
    ggplot2::scale_fill_manual(values = colors) +
    
    # 4. Configure X-axis (time points)
    ggplot2::scale_x_continuous(
      breaks = 1:base::length(base::unique(df$x_pos)),
      labels = base::levels(df$time_point),
      limits = base::c(0.5, base::max(df$x_pos)+0.5)
    ) +
    
    # 5. Configure Y-axis (percentage)
    ggplot2::scale_y_continuous(
      breaks = base::seq(0, 100, 25),
      labels = base::paste0(base::seq(0, 100, 25), "%"),
      limits = base::c(0, 100)
    ) +
    
    # 6. Apply theme styling (preserved original theme)
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # Font settings
      text = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12, color = "black"),
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 13, color = "black"),
      
      # Legend settings
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      legend.key.size = ggplot2::unit(0.8, "cm"),
      
      # Grid settings
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.5),
      
      # Background settings
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Title settings
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 20)),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5,
                                            color = "gray60", margin = ggplot2::margin(b = 20))
    ) +
    
    # 7. Add plot titles and labels (preserved original text)
    ggplot2::labs(
      title = base::paste0("Top 10 over time"),
      subtitle = "Proportion of top 10 integration sites at each time point",
      x = "Time",
      y = "Clonal proportion (%)",
      fill = ""
    )
  
  #################### Optional: Add Value Labels to Stacked Bars ####################
  add_value_labels <- function(p, df) {
    p + ggplot2::geom_text(data = df %>% dplyr::filter(percentage > 5),  # Only show labels for >5% segments
                           ggplot2::aes(x = x_pos,
                                        y = (cum_percentage + cum_percentage_lag) / 2,
                                        label = base::paste0(base::round(percentage, 1), "%")),
                           size = 3, color = "white", fontface = "bold")
  }
  
  # Add value labels (enabled by default, preserved original logic)
  p_with_labels <- add_value_labels(p, df)
  return(p_with_labels)
}
