#' Annotate integration site (IS) data with genomic features
#' 
#' This function adds genomic feature annotations (gene/exon/intron overlap, nearest gene info)
#' to raw integration site data, standardizes chromosome naming, and calculates clone contribution.
#' 
#' @param IS_raw Data frame containing raw IS data with columns: Sample, Chr, Locus, SCount, Strand
#' @return Data frame with annotated genomic features and clone contribution
#' @importFrom GenomicRanges GRanges findOverlaps distanceToNearest
#' @importFrom GenomicFeatures genes exons
#' @export
get_feature=function(IS_raw){
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    stop("Package \"TxDb.Hsapiens.UCSC.hg38.knownGene\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package \"org.Hs.eg.db\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  ## Add 'chr' prefix to Chr column if missing
  if(!base::grepl('chr', IS_raw$Chr)[1]){
    IS_raw$Chr=base::paste0('chr', IS_raw$Chr)
  }
  IS_raw$Chr[IS_raw$Chr=='chr23']='chrX'
  IS_raw$Chr[IS_raw$Chr=='chr24']='chrY'
  IS_raw$Chr[IS_raw$Chr=='chry']='chrY'
  IS_raw$Chr[IS_raw$Chr=='chrx']='chrY'
  
  IS_raw=base::split(IS_raw, f=IS_raw$Sample)
  for(i in base::names(IS_raw)){
    IS_raw[[i]]$Clone_contribution=IS_raw[[i]]$SCount/base::sum(IS_raw[[i]]$SCount)
  }
  IS_raw=base::do.call(base::rbind, IS_raw)
  
  ## Create GRanges object target_gr for target loci
  target_gr <- GenomicRanges::GRanges(
    seqnames = IS_raw$Chr,
    ranges = IRanges::IRanges(start = IS_raw$Locus, end = IS_raw$Locus),  # Single locus
    strand = IS_raw$Strand  # Strand not specified (genes may be on +/-)
  )
  
  ## Load hg38 TxDb annotation
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  ## Get gene/exon info from txdb (suppress messages)
  base::suppressMessages(genes_hg38 <- GenomicFeatures::genes(txdb))
  exons_hg38 <- GenomicFeatures::exons(txdb)
  
  ## 1. Check overlap with gene regions
  gene_overlap <- GenomicRanges::findOverlaps(target_gr, genes_hg38)
  ## Remove duplicate query hits from gene_overlap
  gene_overlap=gene_overlap[!base::duplicated(S4Vectors::queryHits(gene_overlap))]
  
  ## 2. Check overlap with exon regions (filter to gene-overlapping loci)
  exon_overlap <- GenomicRanges::findOverlaps(target_gr, exons_hg38)
  exon_overlap <- exon_overlap[S4Vectors::queryHits(exon_overlap) %in% S4Vectors::queryHits(gene_overlap),]
  
  ## 3. Infer intronic loci (in gene but not in exon)
  IS_raw$in_gene=FALSE
  IS_raw$in_exon=FALSE
  IS_raw$in_gene[S4Vectors::queryHits(gene_overlap)]=TRUE
  IS_raw$in_exon[S4Vectors::queryHits(exon_overlap)]=TRUE
  IS_raw$in_intron=(IS_raw$in_gene==TRUE & IS_raw$in_exon==FALSE)
  
  ## Calculate distance to nearest gene
  nearest_hits=GenomicRanges::distanceToNearest(x = target_gr, subject = genes_hg38)
  ## Get indices/distances of nearest genes
  nearest_gene_idx <-S4Vectors::subjectHits(nearest_hits)
  nearest_distance <- S4Vectors::mcols(nearest_hits)$distance
  
  ## Add nearest gene info to IS_raw
  IS_raw$nearest_gene_idx=nearest_gene_idx
  IS_raw$nearest_distance=nearest_distance
  
  ## Convert genes_hg38 to data frame
  df_genes_hg38=base::data.frame(genes_hg38)
  df_genes_hg38$strand=base::as.character(df_genes_hg38$strand)
  
  ## Add nearest gene ID
  IS_raw$nearest_gene_id=df_genes_hg38$gene_id[IS_raw$nearest_gene_idx]
  
  ## Preserve Strand info and remove original Strand column
  IS_raw$IS_strand=IS_raw$Strand
  IS_raw$Strand=NULL
  
  ## Merge IS_raw with gene info by nearest_gene_id
  IS_raw=base::merge(IS_raw, df_genes_hg38, by.x='nearest_gene_id', by.y='gene_id', all.x=TRUE)
  
  ## Calculate distance to gene start (TSS)
  IS_raw$nearest_distance=base::ifelse(IS_raw$Locus<IS_raw$start,
                                       yes = IS_raw$nearest_distance*-1,
                                       no = IS_raw$nearest_distance)
  
  ## Adjust distance to TSS for genes on '-' strand
  IS_raw$nearest_distance=base::ifelse(IS_raw$strand=='-',
                                       yes = IS_raw$nearest_distance*-1,
                                       no = IS_raw$nearest_distance)
  IS_raw$seqnames=NULL
  
  ## Map Entrez IDs to gene symbols
  IS_raw$nearest_gene_name <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = IS_raw$nearest_gene_id,  # Input: Entrez IDs
    keytype = "ENTREZID",     # Type of input IDs
    column = "SYMBOL",        # Output: gene names (Symbols)
    multiVals = "first"       # Handle multiple matches by taking the first
  )
  
  ### Replace NA gene symbols with Entrez ID
  I=base::is.na(IS_raw$nearest_gene_name)
  IS_raw$nearest_gene_name[I]=IS_raw$nearest_gene_id[I]
  
  ## Remove redundant columns
  IS_raw$nearest_gene_id=NULL
  IS_raw$nearest_gene_idx=NULL
  
  return(IS_raw)
}

