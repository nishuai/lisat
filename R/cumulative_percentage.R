#' Calculate normalized cumulative sum for top N elements of a numeric vector
#' @param x Non-empty numeric vector (integration site ratio data)
#' @return Named vector of cumulative sums for predefined target indices + total sum (all = 1)
#' @import tidyr
#' @export

fit_cum_simple <- function(x) {

  # Validate input: non-empty numeric vector
  if (!base::is.numeric(x) || base::length(x) == 0) {
    base::stop("Input must be a non-empty numeric vector.")
  }
  
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required. Please reinstall the package.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to generate cumulative sum curves. Install with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  if (!requireNamespace("broom", quietly = TRUE)) {
    stop("Package 'broom' is required to generate cumulative sum curves. Install with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  # Sort vector in descending order (highest ratios first)
  x_sorted <- base::sort(x, decreasing = TRUE)
  # Normalize values to sum to 1 (proportion)
  x_norm <- x_sorted / base::sum(x_sorted)
  # Get length of normalized vector
  n <- base::length(x_norm)
  
  # Predefined target indices for cumulative sum calculation
  targets <- base::c(1,2,3, 5,10, 50, 250,1000, 2000)
  # Calculate cumulative sum at target indices (cap at vector length)
  result <- base::cumsum(x_norm)[base::pmin(targets, n)]
  # Add total cumulative sum (all elements = 1)
  result <- base::c(result, "all" = 1)
  # Name result elements with target indices (except "all")
  base::names(result)[-base::length(result)] <- base::as.character(targets)
  
  return(result)
}

#' Plot cumulative curve and perform statistical analysis
#' 
#' @param IS_ratio A numeric vector of integration site ratios (output of fit_cum_simple)
#' @return A list containing the ggplot object, t-test results, and Wilcoxon test result.
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme element_text
#' @importFrom purrr map_dfr map2_dbl
#' @importFrom broom tidy
#' @importFrom stats t.test wilcox.test
#' @importFrom utils combn
#' @export
Cumulative_curve=function(IS_ratio){
  
  index <- NULL
  cum_sum <- NULL
  Group <- NULL
  ### Extract pre-calculated cumulative IS matrix from EPS (internal package object)
  # Ensure EPS is available (it should be in sysdata.rda)
  if (!exists("EPS")) {
      stop("Internal data 'EPS' not found.")
  }
  df_combined=EPS$Cumulative_IS
  
  ### Calculate cumulative sum for input IS_ratio using fit_cum_simple
  guest_data=fit_cum_simple(IS_ratio)
  # Add guest data as new column to combined dataframe
  df_combined=base::cbind(df_combined, guest_data)
  # Duplicate guest data column (for consistency with original code)
  df_combined$guest_data=guest_data
  
  # 2. Reshape data from wide to long format for visualization
  df_long <- df_combined %>%
    tidyr::pivot_longer(
      cols = -index,
      names_to = "vector",
      values_to = "cum_sum"
    )
  
  ### Assign group labels to vectors (Longterm/Cell_Infusion/Sample)
  df_long$Group='Longterm'
  df_long$Group[df_long$vector %in% base::paste0('V', 1:10)]='Cell_Infusion'
  ### Label input guest data as "Sample" group
  df_long$Group[df_long$vector %in% 'guest_data']='Sample'
  
  ## Generate cumulative sum trajectory plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = index, y = cum_sum, color = Group))
  p <- p + ggplot2::geom_line(ggplot2::aes(group = vector), linewidth = 0.8, linetype='dashed')
  p <- p + ggplot2::labs(
    title = "Cumulative Sum Trajectories by Group",
    subtitle = "Comparison between CF and longterm groups",
    x = "Index",
    y = "Cumulative Sum",
    color = "Group"
  )
  p <- p + ggplot2::theme_minimal()
  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 10)
  )
  
  ############# Statistical testing (t-test for each row)
  result_df <- purrr::map_dfr(1:(base::nrow(df_combined)-2), function(row_index) {
    # Extract V1-V20 columns for current row as numeric sample data
    sample_data <- base::as.numeric(df_combined[row_index, base::paste0("V", 1:20)])
    # Extract hypothesized mean (guest_data) for current row
    hypothesized_mean <-df_combined$guest_data[row_index]
    # One-sample t-test: H0 = sample mean equals hypothesized mean
    test_result <- stats::t.test(x = sample_data, mu = hypothesized_mean)
    # Convert test results to tidy dataframe
    broom::tidy(test_result)
  })
  
  # 2. Calculate L2 (Euclidean) distance between curves
  # Step 1: Define L2 distance function
  L2d=function(curve_a, curve_b){
    differences <- curve_a - curve_b
    squared_differences <- differences^2
    sum_of_squares <- base::sum(squared_differences)
    l2_distance <- base::sqrt(sum_of_squares)
    return(l2_distance)
  }
  
  # Test L2 distance calculation (example: column 1 vs column 2)
  L2d(df_combined[,1], df_combined[,2])
  
  # Get columns to compare (V1-V20)
  cols_to_compare <- base::colnames(df_combined)[1:20]
  # Generate all pairwise combinations of columns
  combinations <- utils::combn(cols_to_compare, m = 2, simplify = TRUE)
  # Calculate L2 distance for all column pairs
  distances <- purrr::map2_dbl(
    .x = combinations[1, ],
    .y = combinations[2, ],
    .f = function(col1, col2) {
      # Extract column vectors from dataframe
      vec1 <- df_combined[, col1]
      vec2 <- df_combined[, col2]
      # Calculate L2 distance between two vectors
      return(L2d(vec1, vec2))
    }
  )
  
  # Calculate L2 distance between sample (guest_data) and each V1-V20 column
  Sampe_l2d=c()
  for(i in 1:20){
    Sampe_l2d[i]= L2d(df_combined[,22], df_combined[,i])
  }
  
  # Wilcoxon rank sum test (compare sample L2 distances vs pairwise distances)
  mw_result <- stats::wilcox.test(Sampe_l2d, distances, alternative = "two.sided")
  
  # Return plot + statistical results
  return(base::list(p, result_df, mw_result))
}
