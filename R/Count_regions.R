#' Calculate regional distribution percentages of integration sites (IS)
#' @param IS_raw Data frame of raw integration site data (must have Sample column + regional annotation columns)
#' @param Patient_timepoint Data frame mapping Sample_ID to Time_Point (columns: Sample_ID, Time_Point)
#' @return List of data frames (per sample) with regional IS percentages (Exonic/Intronic/Enhancer etc.)
#' @export
Count_regions <- function(IS_raw, Patient_timepoint) {
  # Initialize list to store regional data per sample
  Region_sample <- list()
  
  ####################################################################
  # Step 1: Split IS data by sample and calculate regional counts
  ####################################################################
  # Split raw IS data into sample-specific subsets
  IS_raw_split <- base::split(IS_raw, f = IS_raw$Sample)
  
  # Iterate over each sample to calculate regional percentages
  for (i in base::names(IS_raw_split)) {
    # Extract sample-specific IS data
    my_data <- IS_raw_split[[i]]
    
    # Define target genomic regions for analysis
    Region_names <- base::c("Exonic", "Intronic", "Integenic", "Enhancer", "Promotor", "SafeHarbor")
    
    # Initialize data frame for regional metrics
    Region_data <- base::data.frame(
      Product = Region_names, 
      Share = NA, 
      Percentage = NA, 
      Label = NA
    )
    
    # Calculate regional counts (proportion of total IS per sample)
    Exonic_count <- base::sum(my_data$in_exon) / base::nrow(my_data)
    Intronic_count <- base::sum(my_data$in_intron) / base::nrow(my_data)
    Integenic_count <- base::sum(!my_data$in_gene) / base::nrow(my_data)
    Enhancer_count <- base::sum(my_data$Enhancer) / base::nrow(my_data)
    Promotor_count <- base::sum(my_data$Promotor) / base::nrow(my_data)
    Safeharbor_count <- base::sum(my_data$Safeharbor) / base::nrow(my_data)
    
    # Assign calculated shares to data frame
    Region_data$Share <- base::c(
      Exonic_count, Intronic_count, Integenic_count, 
      Enhancer_count, Promotor_count, Safeharbor_count
    )
    
    # Convert share to percentage (rounded to 1 decimal place)
    Region_data$Percentage <- base::round(Region_data$Share * 100, 1)
    # Map product names to labels (for visualization compatibility)
    Region_data$Label <- Region_data$Product
    
    # Add sample-specific data to list (temporary Time column with sample name)
    Region_sample[[i]] <- Region_data
    Region_sample[[i]]$Time <- i
  }
  
  ####################################################################
  # Step 2: Map sample names to time points (from Patient_timepoint)
  ####################################################################
  # Iterate over samples to update Time column with actual time points
  for (i in 1:base::length(base::names(Region_sample))) {
    sample_name <- base::names(Region_sample)[i]
    # Check if sample exists in Patient_timepoint mapping
    if (base::sum(Patient_timepoint$Sample_ID == sample_name) == 1) {
      # Update Time column with matched time point
      Region_sample[[i]]$Time <- Patient_timepoint$Time_Point[Patient_timepoint$Sample_ID == sample_name]
    }
  }
  
  # Return list of sample-specific regional distribution data frames
  return(Region_sample)
}