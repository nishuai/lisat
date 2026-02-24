#' Plot chromosome distribution of integration sites (IS)
#' @param IS_raw Data frame containing raw integration site data (must have Chr column)
#' @param ref_version Reference version for simulation (options: 'random' or 'LV', default = 'random')
#' @return ggplot object of chromosome distribution (percentage of IS per chromosome)
#' @export
chr_distribution <- function(IS_raw, ref_version = 'random') {
  # Convert chromosome column to factor with standard human chromosome levels
  Freq <- NULL
  variable <- NULL
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to use this function. Please install it with:\ninstall.packages('ggplot2')", call. = FALSE)
  }
  chrs <- IS_raw$Chr
  chrs <- base::as.character(chrs)
  # Set factor levels to standard human chromosomes (1-22, X, Y, M)
  chrs <- base::factor(chrs, levels = base::paste0('chr', c(base::as.character(1:22), 'X','Y','M')))
  # Count IS frequency per chromosome
  chrs <- base::data.frame(base::table(chrs))
  
  # Normalize frequency to percentage (per chromosome)
  chrs$Freq <- (chrs$Freq / base::sum(chrs$Freq)) * 100
  
  # Assign simulation frequency based on reference version
  if (base::tolower(ref_version) == 'random') {
    simulation <- chrs
    # Predefined simulation frequencies for 'random' reference
    simulation$Freq <- c(8.030, 7.804, 6.389, 6.145, 5.884, 5.536, 5.196, 4.722,
                         4.463, 4.304, 4.367, 4.345, 3.728, 3.481, 3.299, 2.934,
                         2.694, 2.618, 1.891, 2.082, 1.515, 1.650, 5.066, 1.858,
                         0.00053)
  }
  
  if (base::tolower(ref_version) == 'LV') {
    simulation <- chrs
    # Predefined simulation frequencies for 'LV' reference
    simulation$Freq <- c(8.098, 8.203, 6.959, 5.968, 5.736, 5.884, 5.799, 5.377,
                         4.471, 4.513, 4.660, 5.124, 3.458, 2.973, 2.847, 2.720,
                         2.509, 2.741, 2.024, 2.067, 1.476, 0.970, 4.934, 0.485,
                         0.00053)
  }
  
  # Prepare data for plotting (merge experimental and simulation data)
  my_df <- base::rbind(chrs, simulation)
  my_df$variable <- base::rep(c('Exp', 'in_silico'), each = 25)
  my_df$variable <- base::factor(my_df$variable, levels = c('Exp','in_silico'))
  
  # Generate chromosome distribution plot
  p <- ggplot2::ggplot(data = my_df, ggplot2::aes(x = chrs, y = Freq, fill = variable)) +
    ggplot2::geom_bar(stat = 'identity', position = 'dodge', width = 0.66) +
    ggplot2::scale_fill_manual(values = c('steelblue','grey70')) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::labs(fill = '') +
    ggplot2::xlab('Chromosome') +
    ggplot2::ylab('IS in each chromosome [%]')
  
  return(p)
}