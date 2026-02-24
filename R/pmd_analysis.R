#' Calculate PMD (Proportional Modular Diversity) for integration site data
#' 
#' This function computes UIS count, top clone contribution percentage, and PMD metrics (Richness/Eveness/PMD)
#' for integration site data, and maps samples to patient time points.
#' 
#' @param IS_raw Data frame containing integration site data (columns: Sample, Clone_contribution required)
#' @param Patient_timepoint Data frame mapping Sample_ID to Time_Point (columns: Sample_ID, Time_Point required)
#' @return Data frame with PMD metrics (UIS, TOP_P, Richness, Eveness, PMD, Sample, Time)
#' @export
pmd_analysis=function(IS_raw, Patient_timepoint){
  
  IS_raw=base::split(IS_raw, f=IS_raw$Sample)
  UIS_count=base::sapply(IS_raw, function(x) base::nrow(x))
  UIS_top1=base::sapply(IS_raw, function(x) base::max(x$Clone_contribution))*100
  
  my_df=base::data.frame(UIS=UIS_count, TOP_P=UIS_top1)
  my_df$UIS=base::as.numeric(my_df$UIS)
  my_df$Richness=base::log2(my_df$UIS)
  my_df$Eveness=base::log2(1/(my_df$TOP_P/100))
  my_df$PMD=base::log2(my_df$Eveness/(my_df$Richness-my_df$Eveness))
  my_df$Sample=base::rownames(my_df)
  
  PMD_data=my_df
  PMD_data$Time=NA
  
  for(i in 1:base::nrow(PMD_data)){
    if(base::sum(Patient_timepoint$Sample_ID==PMD_data$Sample[i])==1){
      PMD_data$Time[i]=Patient_timepoint$Time_Point[Patient_timepoint$Sample_ID==PMD_data$Sample[i]]
    }
  }
  
  if(base::any(base::is.na(PMD_data$Time))){
    PMD_data$Time=base::sapply(base::strsplit(PMD_data$Sample, '-'),function(x) x[2])
  }
  
  natural_sort <- function(x) {
    numeric_part <- base::as.numeric(base::gsub("[^0-9.]", "", x))
    unit_part <- base::gsub("[0-9.]", "", x)
    base::order(numeric_part, unit_part)
  }
  
  PMD_data$Time=base::factor(PMD_data$Time, levels=PMD_data$Time[natural_sort(PMD_data$Time)])
  PMD_data=PMD_data[natural_sort(PMD_data$Time),]
  
  return(PMD_data)
}
