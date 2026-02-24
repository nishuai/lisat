#' Visualize and analyze network of common integration sites (CIS)
#' @param IS_raw Data frame containing integration site data (must have Locus, Chr, nearest_gene_name columns)
#' @param connect_distance Numeric threshold for connecting IS (default = 50000 bp)
#' @return Data frame with top 10 CIS network metrics (Chr, Locus, Gene, Total_dots, etc.)
#' @importFrom utils head data
#' @importFrom stats dist kmeans
#' @export
CIS <- function(IS_raw, connect_distance = 50000) {
  
  # Print function status message
  base::cat('CIS_network calculates interactions between integration sites in a network view\n\n')
  
  ####################################################################
  # Step 1: Validate and preprocess input data
  ####################################################################
  # Extract core columns for CIS analysis
  cis_data <- IS_raw[, c('Locus','Chr','nearest_gene_name')]
  # Remove rows with missing chromosome information
  cis_data <- cis_data[!base::is.na(cis_data$Chr), ]
  
  # Validate input column requirements
  if (!base::ncol(cis_data) == 3) {
    base::stop('Input data must contain exactly 3 columns: "Locus", "Chr", "nearest_gene_name"')
  }
  
  # Remove duplicate records
  cis_data <- cis_data[!(base::duplicated(cis_data)), ]
  # Split data by chromosome for per-chromosome analysis
  cis_data_split <- base::split(cis_data, f = cis_data$Chr)
  
  ####################################################################
  # Step 2: Build network edges between integration sites
  ####################################################################
  # Initialize edge data frame
  edges <- base::data.frame(from = NA, to = NA, width = NA)
  
  # Iterate over each chromosome to identify edges
  for (k in 1:base::length(cis_data_split)) {
    my_df <- cis_data_split[[k]]
    # Create unique ID for each integration site (gene + index + chr + locus)
    my_df$nearest_gene_name <- base::paste0(
      my_df$nearest_gene_name, '___', 1:base::nrow(my_df),
      '__', my_df$Chr, '__', my_df$Locus
    )
    # Sort by locus position
    my_df <- my_df[base::order(my_df$Locus), ]
    
    # Calculate distance matrix and filter by connection threshold
    dist_matrix <- base::as.matrix(stats::dist(my_df$Locus))
    dist_matrix <- dist_matrix < connect_distance
    # Remove upper triangle (avoid duplicate edges)
    for(row in 1:base::nrow(dist_matrix)){
      dist_matrix[row, 1:row] <- FALSE
    }
    
    # Extract edge pairs from distance matrix
    edge_left <- base::rep(my_df$nearest_gene_name, base::rowSums(dist_matrix))
    edge_right <- c()
    for (i in 1:base::nrow(dist_matrix)){
      edge_right <- base::c(edge_right, my_df$nearest_gene_name[dist_matrix[i,]])
    }
    
    # Add valid edges to edge data frame
    if(base::length(edge_left) > 0){
      edge <- base::data.frame(from = edge_left, to = edge_right, width = 1)
      edges <- base::rbind(edges, edge)
    }
  }
  
  # Remove initial dummy row from edge data frame
  edges <- edges[-1, ]
  
  # Exit early if no valid edges (save empty result)
  if(base::nrow(edges) <= 1){
    base::stop('No valid edges found')
  }
  
  ####################################################################
  # Step 3: Filter to top genes and refine network
  ####################################################################
  # Extract gene names from edge IDs
  all_nodes <- base::c(edges$from, edges$to)
  all_nodes_gene <- base::sapply(base::strsplit(all_nodes, '___'), function(x) x[1])
  # Get top 30 genes by connection frequency
  top_50 <- head(base::sort(base::table(all_nodes_gene), decreasing = TRUE), 30)
  
  # Filter edges to include only top 30 genes
  edges_from_gene <- base::sapply(base::strsplit(edges$from, '___'), function(x) x[1])
  edges_to_gene <- base::sapply(base::strsplit(edges$to, '___'), function(x) x[1])
  edges <- edges[edges_from_gene %in% base::names(top_50) |
                   edges_to_gene %in% base::names(top_50) , ]
  
  # Identify representative nodes for top genes
  nodes_rep <- base::data.frame(id = base::c(edges$from, edges$to), label = NA)
  nodes_rep$label <- base::sapply(base::strsplit(nodes_rep$id, '___'), function(x) x[1])
  
  top_net <- list()
  for (i in 1:base::length(top_50)){
    # Get most representative node for each top gene
    top_reps <- nodes_rep$id[nodes_rep$label == base::names(top_50)[i]]
    top_net[[i]] <- base::names(base::sort(base::table(top_reps), decreasing = TRUE)[1])
  }
  
  # Assign edges to top gene networks (expand connections iteratively)
  edges$belong <- ''
  for (k in 1:base::length(top_net)){
    # Iterate to expand network connections (20 iterations for indirect connections)
    for (do in 1:20){
      for (i in 1:base::nrow(edges)){
        if (edges$from[i] %in% top_net[[k]] | edges$to[i] %in% top_net[[k]]){
          edges$belong[i] <- base::names(top_50)[k]
          top_net[[k]] <- base::unique(base::c(top_net[[k]], base::c(edges$from[i], edges$to[i])))
        }
      }
    }
  }
  
  # Filter to top 10 largest networks
  edges <- edges[edges$belong != '', ]
  top_10 <- base::names(base::sort(base::table(edges$belong), decreasing = TRUE))[1:base::min(10, base::length(base::unique(edges$belong)))]
  edges <- edges[edges$belong %in% top_10, ]
  
  # Prepare node labels (limit to 1 label per gene to avoid overlap)
  nodes <- base::data.frame(id = base::unique(base::c(edges$from, edges$to)), label = NA)
  nodes$label <- base::sapply(base::strsplit(nodes$id, '___'), function(x) x[1])
  
  # Randomly select 1 representative label per gene (set seed for reproducibility)
  base::set.seed(100)
  nodes <- nodes[base::order(nodes$label), ]
  my_index <- base::split(nodes$label, f = nodes$label)
  
  index_length <- c()
  for (i in 1:base::length(my_index)){
    index_length <- base::c(index_length, base::length(my_index[[i]]))
    my_index[[i]] <- 1:base::length(my_index[[i]])
    if(base::length(my_index[[i]]) > 1){
      my_index[[i]] <- base::sample(my_index[[i]], 1)
    }
  }
  
  # Update label visibility (hide duplicate labels)
  index_length <- base::cumsum(index_length)
  for (i in 2:base::length(my_index)){
    my_index[[i]] <- my_index[[i]] + index_length[i-1]
  }
  my_index <- base::unlist(my_index)
  nodes$keep <- FALSE
  nodes$keep[my_index] <- TRUE
  nodes$label[base::which(nodes$keep == FALSE)] <- ''
  nodes$keep <- NULL
  
  # Build igraph network object (undirected)
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  ####################################################################
  # Step 4: Calculate CIS network metrics for top 10 networks
  ####################################################################
  tops <- list()
  for (i in 1:base::length(top_10)){
    gene_base <- top_10[i]
    # Extract chromosome information
    chr <- base::unique(base::sapply(base::strsplit(edges$from[edges$belong == gene_base], '__'), function(x) x[3]))
    # Extract and convert locus positions
    all_locus <- base::sapply(base::strsplit(edges$from[edges$belong == gene_base], '__'), function(x) x[4])
    all_locus <- base::as.numeric(all_locus)
    # Identify central locus (most frequent)
    Locus <- base::names(base::sort(base::table(all_locus), decreasing = TRUE))[1]
    Locus <- base::as.numeric(Locus)
    
    # Calculate network metrics
    Central_connectivity <- base::as.numeric(base::table(all_locus)[1])
    Dimention <- base::abs(base::range(all_locus)[1] - base::range(all_locus)[2])
    order <- base::sum(edges$belong == gene_base) + 1
    
    # Extract connected genes
    genes_from <- base::sapply(base::strsplit(edges$from[edges$belong == gene_base], '__'), function(x) x[1])
    genes_to <- base::sapply(base::strsplit(edges$to[edges$belong == gene_base], '__'), function(x) x[1])
    gene_network <- base::unique(base::c(genes_from, genes_to))
    gene_network <- base::c(gene_base, gene_network)
    gene_network <- base::paste(base::unique(gene_network), collapse = '; ')
    
    # Calculate distance metrics
    from_pos <- base::as.numeric(base::sapply(base::strsplit(edges$from[edges$belong == gene_base], '__'), function(x) x[4]))
    to_pos <- base::as.numeric(base::sapply(base::strsplit(edges$to[edges$belong == gene_base], '__'), function(x) x[4]))
    total_dots <- base::length(base::unique(base::c(from_pos, to_pos)))
    mean_dist <- base::mean(base::abs(from_pos - to_pos))
    
    # Calculate splitness (kmeans clustering metric)
    if(base::length(from_pos - to_pos) > 2){
      splitness <- stats::kmeans(base::abs(from_pos - to_pos), centers = 2)$betweenss /
        stats::kmeans(base::abs(from_pos - to_pos), centers = 2)$totss
    } else {
      splitness <- 1
    }
    
    # Calculate cluster index
    cluster_index <- base::round(order / (total_dots * (total_dots - 1) / 2), 2)
    
    # Store metrics for current top gene
    tops[[i]] <- base::c(chr, Locus, gene_base, total_dots, Central_connectivity,
                         Dimention, order, gene_network, mean_dist, splitness,
                         cluster_index)
  }
  
  # Convert metrics list to data frame and format
  tops <- base::data.frame(base::do.call(rbind, tops))
  tops$TOP <- base::paste('Top', 1:base::min(base::nrow(tops), 10))
  base::names(tops) <- c('Chr','Locus','Gene','Total_dots','Central_connectivity',
                         'Dimention','Order','Gene_network','Mean_dist','Splitness',
                         'Cluster_index','Top10 CIS')
  
  # Convert numeric columns to appropriate type
  tops$Locus <- base::as.numeric(tops$Locus)
  tops$Dimention <- base::as.numeric(tops$Dimention)
  tops$Order <- base::as.numeric(tops$Order)
  # Clean gene names (replace underscores with spaces)
  tops$Gene <- base::gsub('_', ' ', tops$Gene)
  
  return(tops)
}
