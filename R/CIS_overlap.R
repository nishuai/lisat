#' Generate colored GT table for CIS overlap across samples/timepoints
#' @param CIS_data Data frame of CIS metrics (must have Chr and Locus columns)
#' @param IS_raw Data frame of raw integration site data (must have Sample, Chr, Locus columns)
#' @param Timelevels Optional vector of sample/timepoint levels for ordered display (default = NULL)
#' @return gt table object with colored CIS overlap status (TRUE/FALSE)
#' @importFrom stats na.omit
#' @export
CIS_overlap <- function(CIS_data, IS_raw, Timelevels = NULL) {
  
  ####################################################################
  # Step 1: Define color palette for table formatting
  ####################################################################
  # Base color palette + RColorBrewer palettes (expanded for flexibility)
  colors <- base::c(
    "grey60","#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4",
    "#FECA57", "#FF9FF3", "#54A0FF", "#5F27CD"
  )
  # Add additional brewer palettes
  colors <- base::c(colors, RColorBrewer::brewer.pal('Set2', n = 8))
  colors <- base::c(colors, RColorBrewer::brewer.pal('Set1', n = 8))
  colors <- base::c(colors, RColorBrewer::brewer.pal('RdBu', n = 9))
  # Remove duplicate/unwanted colors
  colors <- colors[-(30:31)]
  
  ####################################################################
  # Step 2: Preprocess CIS and sample data
  ####################################################################
  # Extract core CIS metrics (first 8 columns)
  CIS_overtime <- CIS_data[, 1:8]
  # Split integration site data by sample/timepoint
  my_sample <- base::split(IS_raw, f = IS_raw$Sample)
  
  ####################################################################
  # Step 3: Define natural sort function for sample/timepoint names
  ####################################################################
  natural_sort <- function(x) {
    # Extract numeric and non-numeric parts for natural sorting
    numeric_part <- base::as.numeric(base::gsub("[^0-9.]", "", x))
    unit_part <- base::gsub("[0-9.]", "", x)
    # Return sorted index
    base::order(numeric_part, unit_part)
  }
  
  # Sort samples (natural sort if Timelevels not provided)
  if (base::is.null(Timelevels)) {
    my_sample <- my_sample[natural_sort(base::names(my_sample))]
  } else {
    my_sample <- my_sample[Timelevels]
  }
  
  # Add empty columns for each sample (to store overlap status)
  for (i in base::names(my_sample)) {
    CIS_overtime$dummy <- NA
    base::names(CIS_overtime)[base::ncol(CIS_overtime)] <- i
  }
  
  ####################################################################
  # Step 4: Check CIS overlap in each sample/timepoint
  ####################################################################
  # Iterate over each CIS entry
  for (i in 1:base::nrow(CIS_overtime)) {
    # Get query chromosome and locus for current CIS
    chr_query <- CIS_overtime$Chr[i]
    loci_query <- CIS_overtime$Locus[i]
    
    # Check overlap in each sample (50kb threshold)
    for (k in base::names(my_sample)) {
      # Calculate locus distance and check if <50kb (CIS overlap)
      CIS_exist <- base::sum(
        my_sample[[k]]$Chr == chr_query &
          base::abs(base::as.numeric(loci_query) - base::as.numeric(my_sample[[k]]$Locus)) < 50000
      ) > 0
      # Assign overlap status (TRUE/FALSE) to sample column
      CIS_overtime[i, k] <- CIS_exist
    }
  }
  
  ####################################################################
  # Step 5: Generate colored GT table
  ####################################################################
  # Create base GT table object
  gt_table <- gt::gt(CIS_overtime)
  
  # Apply color coding to sample columns (TRUE = colored, FALSE = grey)
  for (i in base::names(my_sample)) {
    gt_table <- gt::data_color(
      gt_table,
      columns = i,
      fn = scales::col_factor(
        palette = c("FALSE" = colors[1], "TRUE" = colors[2]),
        levels = c("FALSE","TRUE")
      )
    )
  }
  
  # Return formatted GT table (for display in RStudio/export)
  return(gt_table)
}