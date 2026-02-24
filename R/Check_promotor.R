#' Check if integration sites (IS) are located in promoter regions
#' @param IS_raw Data frame containing raw integration site data (must have Chr and Locus columns)
#' @return Data frame with an added Promotor column (TRUE = located in promoter, FALSE = not located in promoter)
#' @export
Promotor_check <- function(IS_raw) {
  # Load internal promoter reference data (EPS is an internal package object)
  promotor <- EPS[['promotor']]
  # Split promoter region data by chromosome
  promotor_chr <- base::split(promotor, f = promotor$V1)
  
  # Unify chromosome naming (add chr prefix to avoid matching failure)
  if (!base::grepl('chr', IS_raw$Chr[1])) {
    IS_raw$Chr <- base::paste0('chr', IS_raw$Chr)
  }
  # Split integration site data by chromosome
  dd1 <- base::split(IS_raw, f = IS_raw$Chr)
  promotor_loci <- c()  # Store loci located in promoter regions
  
  # Match promoter regions chromosome by chromosome
  for (i in base::names(dd1)) {
    # Merge loci with promoter start/end positions for sorted comparison
    my_df <- base::data.frame(
      loci = base::c(dd1[[i]]$Locus, promotor_chr[[i]]$V2, promotor_chr[[i]]$V3),
      name = base::c(
        base::rep('loci', base::length(dd1[[i]]$Locus)),
        base::rep('start', base::length(promotor_chr[[i]]$V2)),
        base::rep('end', base::length(promotor_chr[[i]]$V3))
      )
    )
    # Sort by locus position to ensure correct comparison order
    my_df <- my_df[base::order(my_df$loci), ]
    in_promotor <- FALSE  # Flag to mark if entering promoter region
    
    # Iterate through loci to judge if in promoter interval
    for (k in 1:base::nrow(my_df)) {
      if (my_df$name[k] == 'start') in_promotor <- TRUE  # Enter promoter start position
      if (my_df$name[k] == 'end') in_promotor <- FALSE    # Leave promoter end position
      # Record integration loci located in promoter intervals
      if (in_promotor == TRUE & my_df$name[k] == 'loci') {
        promotor_loci <- base::c(promotor_loci, base::paste(i, my_df$loci[k]))
      }
    }
  }
  
  # Label if integration sites are in promoter regions
  if (base::length(promotor_loci) > 0) {
    my_sheet_I <- c()
    # Match row indices corresponding to the loci
    for (i in 1:base::length(promotor_loci)) {
      chrs <- base::strsplit(promotor_loci[i], ' ')[[1]][1]
      locis <- base::strsplit(promotor_loci[i], ' ')[[1]][2]
      my_I <- base::which(IS_raw$Chr == chrs & IS_raw$Locus == locis)
      my_sheet_I <- base::c(my_sheet_I, my_I)
    }
    # Add Promotor column and label results
    IS_raw$Promotor <- FALSE
    IS_raw$Promotor[my_sheet_I] <- TRUE
    return(IS_raw)
  } else {
    # Set Promotor column to FALSE when no matching loci
    IS_raw$Promotor <- FALSE
    return(IS_raw)
  }
}