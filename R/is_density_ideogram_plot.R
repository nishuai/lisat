#' Plot chromosome ideogram with integration site annotations
#' 
#' This function generates a chromosome ideogram plot showing the density and position
#' of integration sites (IS) using the RIdeogram package.
#' 
#' @param IS_raw Data frame containing integration site data (columns: Chr, Locus required)
#' @param output_dir Character, path to output directory for the PDF plot
#' @importFrom utils head data
#' @return None (generates a PDF file in output_dir)
#' @export
ideogram_plot=function(IS_raw, output_dir){
  if (!requireNamespace("RIdeogram", quietly = TRUE)) {
    stop("Package 'RIdeogram' is required. Please reinstall the package.", call. = FALSE)
  }
  if (!requireNamespace("Cairo", quietly = TRUE)) {
    stop("Package 'Cairo' is required for ideogram plotting. Please install it.")
  }
  
  human_karyotype <- NULL
  gene_density <- NULL
  Random_RNAs_500 <- NULL
  ####read built-in data from RIdeogram
  data(human_karyotype, package="RIdeogram", envir = environment())
  data(gene_density, package="RIdeogram", envir = environment())
  data(Random_RNAs_500, package="RIdeogram", envir = environment())
  
  ######process input data
  mydf=IS_raw
  mymarks=data.frame(matrix(NA, nrow = base::nrow(mydf), ncol=base::ncol(Random_RNAs_500)))
  base::names(mymarks)=base::names(Random_RNAs_500)
  mymarks$Type='Insertion'
  mymarks$Shape='triangle'
  mymarks$Chr=mydf$Chr
  mymarks$Chr=base::gsub('chr','', mymarks$Chr)  #
  mymarks$Start=mydf$Locus
  mymarks$End=mymarks$Start+100
  
  is_density=gene_density
  is_density$Value=0
  
  mymarks=mymarks[!base::is.na(mymarks$Start),]
  
  for (i in 1:base::nrow(mymarks)){
    I1=mymarks$Chr[i]==is_density$Chr
    I2=mymarks$Start[i]>is_density$Start
    I3=mymarks$Start[i]<is_density$End
    final_loci=base::which(I1 & I2 & I3)
    if (base::length(final_loci) > 0) { 
      is_density$Value[final_loci]=is_density$Value[final_loci]+1
    }
  }
  
  if (!base::dir.exists(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE)  
  }
  
  svg_file <- file.path(output_dir, "chromosome.svg")
  
  Cairo::CairoSVG(svg_file, 
                  width = 10, height = 6)  #  Cairo SVG device
  old_device <- getOption("device") 
  options(device = Cairo::CairoSVG)
  
  RIdeogram::ideogram(karyotype = human_karyotype,
                      overlaid = is_density,
                      label_type = "marker",
                      width = 190,
                      colorset1 = base::c("#fffefe", "#ef284e", "#230108"))
  
  suppressWarnings(
  RIdeogram::convertSVG(svg = svg_file,
                        file=base::file.path(output_dir, 'Ideogram.pdf'), 
                        device = "pdf")
  )

}
