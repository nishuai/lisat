#' Check if integration sites (IS) are located in safe harbor regions
#' @param IS_raw Data frame containing raw integration site data (must have Chr and Locus columns)
#' @return Data frame with an added Safeharbor column (TRUE = located in safe harbor, FALSE = not located in safe harbor)
#' @export
Safeharbor_check <- function(IS_raw) {
  # Load internal safe harbor reference data (EPS is an internal package object)
  safeharbor <- EPS[['safeharbor']]
  # Split safe harbor region data by chromosome
  safeharbor_chr <- base::split(safeharbor, f = safeharbor$V1)
  
  # Unify chromosome naming (add chr prefix to avoid matching failure)
  if (!base::grepl('chr', IS_raw$Chr[1])) {
    IS_raw$Chr <- base::paste0('chr', IS_raw$Chr)
  }
  # Split integration site data by chromosome
  dd1 <- base::split(IS_raw, f = IS_raw$Chr)
  safeharbor_loci <- c()  # Store loci located in safe harbor regions
  
  # Match safe harbor regions chromosome by chromosome
  for (i in base::names(dd1)) {
    # Merge loci with safe harbor start/end positions for sorted comparison
    my_df <- base::data.frame(
      loci = base::c(dd1[[i]]$Locus, safeharbor_chr[[i]]$V2, safeharbor_chr[[i]]$V3),
      name = base::c(
        base::rep('loci', base::length(dd1[[i]]$Locus)),
        base::rep('start', base::length(safeharbor_chr[[i]]$V2)),
        base::rep('end', base::length(safeharbor_chr[[i]]$V3))
      )
    )
    # Sort by locus position to ensure correct comparison order
    my_df <- my_df[base::order(my_df$loci), ]
    in_safeharbor <- FALSE  # Flag to mark if entering safe harbor region
    
    # Iterate through loci to judge if in safe harbor interval
    for (k in 1:base::nrow(my_df)) {
      if (my_df$name[k] == 'start') in_safeharbor <- TRUE  # Enter safe harbor start position
      if (my_df$name[k] == 'end') in_safeharbor <- FALSE    # Leave safe harbor end position
      # Record integration loci located in safe harbor intervals
      if (in_safeharbor == TRUE & my_df$name[k] == 'loci') {
        safeharbor_loci <- base::c(safeharbor_loci, base::paste(i, my_df$loci[k]))
      }
    }
  }
  
  # Label if integration sites are in safe harbor regions
  if (base::length(safeharbor_loci) > 0) {
    my_sheet_I <- c()
    # Match row indices corresponding to the loci
    for (i in 1:base::length(safeharbor_loci)) {
      chrs <- base::strsplit(safeharbor_loci[i], ' ')[[1]][1]
      locis <- base::strsplit(safeharbor_loci[i], ' ')[[1]][2]
      my_I <- base::which(IS_raw$Chr == chrs & IS_raw$Locus == locis)
      my_sheet_I <- base::c(my_sheet_I, my_I)
    }
    # Add Safeharbor column and label results
    IS_raw$Safeharbor <- FALSE
    IS_raw$Safeharbor[my_sheet_I] <- TRUE
    return(IS_raw)
  } else {
    # Set Safeharbor column to FALSE when no matching loci
    IS_raw$Safeharbor <- FALSE
    return(IS_raw)
  }
}

