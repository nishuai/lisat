#' Check if integration sites (IS) are located in enhancer regions
#' @param IS_raw Data frame containing raw integration site data (must have Chr and Locus columns)
#' @return Data frame with an added Enhancer column (TRUE = located in enhancer, FALSE = not located in enhancer)
#' @export
Enhancer_check <- function(IS_raw) {
  # Load internal enhancer reference data (EPS is an internal package object)
  enhancer <- EPS[['enhancer']]
  # Split enhancer region data by chromosome
  enhancer_chr <- base::split(enhancer, f = enhancer$V1)
  
  # Unify chromosome naming (add chr prefix to avoid matching failure)
  if (!base::grepl('chr', IS_raw$Chr[1])) {
    IS_raw$Chr <- base::paste0('chr', IS_raw$Chr)
  }
  # Split integration site data by chromosome
  dd1 <- base::split(IS_raw, f = IS_raw$Chr)
  enhancer_loci <- c()  # Store loci located in enhancer regions
  
  # Match enhancer regions chromosome by chromosome
  for (i in base::names(dd1)) {
    # Merge loci with enhancer start/end positions for sorted comparison
    my_df <- base::data.frame(
      loci = base::c(dd1[[i]]$Locus, enhancer_chr[[i]]$V2, enhancer_chr[[i]]$V3),
      name = base::c(
        base::rep('loci', base::length(dd1[[i]]$Locus)),
        base::rep('start', base::length(enhancer_chr[[i]]$V2)),
        base::rep('end', base::length(enhancer_chr[[i]]$V3))
      )
    )
    # Sort by locus position to ensure correct comparison order
    my_df <- my_df[base::order(my_df$loci), ]
    in_enhancer <- FALSE  # Flag to mark if entering enhancer region
    
    # Iterate through loci to judge if in enhancer interval
    for (k in 1:base::nrow(my_df)) {
      if (my_df$name[k] == 'start') in_enhancer <- TRUE  # Enter enhancer start position
      if (my_df$name[k] == 'end') in_enhancer <- FALSE    # Leave enhancer end position
      # Record integration loci located in enhancer intervals
      if (in_enhancer == TRUE & my_df$name[k] == 'loci') {
        enhancer_loci <- base::c(enhancer_loci, base::paste(i, my_df$loci[k]))
      }
    }
  }
  
  # Label if integration sites are in enhancer regions
  if (base::length(enhancer_loci) > 0) {
    my_sheet_I <- c()
    # Match row indices corresponding to the loci
    for (i in 1:base::length(enhancer_loci)) {
      chrs <- base::strsplit(enhancer_loci[i], ' ')[[1]][1]
      locis <- base::strsplit(enhancer_loci[i], ' ')[[1]][2]
      my_I <- base::which(IS_raw$Chr == chrs & IS_raw$Locus == locis)
      my_sheet_I <- base::c(my_sheet_I, my_I)
    }
    # Add Enhancer column and label results
    IS_raw$Enhancer <- FALSE
    IS_raw$Enhancer[my_sheet_I] <- TRUE
    return(IS_raw)
  } else {
    # Fix: Corrected Promotor column assignment error to Enhancer
    IS_raw$Enhancer <- FALSE
    return(IS_raw)
  }
}
