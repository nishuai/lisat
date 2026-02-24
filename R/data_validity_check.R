#' Validate and standardize integration site (IS) raw data frame
#' @param IS_raw Data frame containing IS data (expected columns: Sample, SCount, Chr, Locus)
#' @return List with validation results: 
#'   - valid (logical): TRUE if data passes validation, FALSE otherwise
#'   - errors (character): Validation messages/errors
#'   - converted_data (data.frame): Original/cleaned data with numeric conversions (if applicable)
#' @export
validate_IS_raw <- function(IS_raw) {
  # Define expected core columns for IS data validation
  expected_cols <- base::c('Sample', 'SCount', 'Chr', 'Locus')
  
  # Initialize validation result list (default to valid = TRUE)
  result <- base::list(
    valid = TRUE,
    errors = base::character(0),
    converted_data = IS_raw  # Store original data (will update with conversions)
  )
  
  # 1. Validate input is a data frame
  if (!base::is.data.frame(IS_raw)) {
    result$valid <- FALSE
    result$errors <- base::c(result$errors, "Input must be a data frame")
    return(result)
  }
  
  # 2. Check for missing required columns
  missing_cols <- base::setdiff(expected_cols, base::colnames(IS_raw))
  if (base::length(missing_cols) > 0) {
    result$valid <- FALSE
    result$errors <- base::c(result$errors,
                             base::paste("Missing required columns:", base::paste(missing_cols, collapse = ", ")))
    return(result)
  }
  
  # 3. Flag extra columns (informational only, not an error)
  extra_cols <- base::setdiff(base::colnames(IS_raw), expected_cols)
  if (base::length(extra_cols) > 0) {
    result$errors <- base::c(result$errors,
                             base::paste("Extra columns present:", base::paste(extra_cols, collapse = ", ")))
  }
  
  # 4. Validate Sample column type (must be character)
  if (!base::is.character(IS_raw$Sample)) {
    result$valid <- FALSE
    result$errors <- base::c(result$errors, "Sample column must be character type")
  }
  
  # 5. Validate and convert strict numeric columns (SCount, Locus)
  strict_numeric_cols <- base::c("SCount", "Locus")
  for (col in strict_numeric_cols) {
    if (base::is.character(IS_raw[[col]])) {
      # Attempt to convert character column to numeric (suppress warning for non-numeric values)
      converted <- base::suppressWarnings(base::as.numeric(IS_raw[[col]]))
      
      # Check if conversion introduced NAs (indicates invalid numeric strings)
      if (base::any(base::is.na(converted) & !base::is.na(IS_raw[[col]]))) {
        result$valid <- FALSE
        result$errors <- base::c(result$errors,
                                 base::paste(col, "column contains non-numeric values that cannot be converted"))
      } else {
        # Successful conversion - update converted_data with numeric values
        result$converted_data[[col]] <- converted
        result$errors <- base::c(result$errors,
                                 base::paste(col, "column was converted from character to numeric"))
      }
    } else if (!base::is.numeric(IS_raw[[col]])) {
      # Column is neither character nor numeric - unrecoverable error
      result$valid <- FALSE
      result$errors <- base::c(result$errors,
                               base::paste(col, "column must be numeric or character with numeric values"))
    }
  }
  
  # 6. Validate Chr column (allow numeric values + X/M/Y (case-insensitive))
  if ("Chr" %in% base::colnames(IS_raw)) {
    chr_col <- IS_raw$Chr
    
    # Convert Chr column to character for consistent validation (preserves X/M/Y)
    if (!base::is.character(chr_col)) {
      chr_col <- base::as.character(chr_col)
    }
    
    # Define valid non-numeric chromosome values (case-insensitive)
    valid_chr_values <- base::c("X", "M", "Y")
    
    # Check each value in Chr column for validity
    for (i in base::seq_along(chr_col)) {
      val <- chr_col[i]
      if (!base::is.na(val)) {
        # Check if value is a valid letter (X/M/Y)
        if (!base::toupper(val) %in% valid_chr_values) {
          # Try to convert to numeric (valid for autosomes 1-22)
          num_val <- base::suppressWarnings(base::as.numeric(val))
          if (base::is.na(num_val)) {
            # Invalid value (not number/X/M/Y)
            result$valid <- FALSE
            result$errors <- base::c(result$errors,
                                     base::paste("Chr column contains invalid value '", val, "' at position", i,
                                                 "- only numbers, X, M, Y are allowed"))
          }
        }
      }
    }
    
    # Preserve validated Chr values (including X/M/Y) in converted_data
    result$converted_data$Chr <- chr_col
  }
  
  # 7. Stop execution if validation failed (return errors)
  if (!result$valid) {
    # Throw error with all validation messages
    base::stop(result$errors)
  }
  
  # Return full validation result (valid + errors + converted data)
  return(result)
}
