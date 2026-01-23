#' Generate Pseudo-Replicates from a Seurat Object by Bootstrap or Splitting
#'
#' Creates pseudo-replicates from single-cell RNA-seq data stored in a Seurat object,
#' using one of two statistical approaches:
#'
#' 1. **Bootstrap sampling ("bootstrap" method):**
#'    This method performs sampling **with replacement** of cells within each sample group.
#'    Each pseudo-replicate contains the same number of cells as the original sample,
#'    but some cells may be repeated and others omitted due to resampling.
#'    Bootstrapping approximates the sampling distribution of the data, allowing
#'    assessment of variability and robustness of downstream analyses.
#'    Because it creates a new count matrix by resampling cells, this method returns
#'    a new Seurat object with freshly sampled counts and metadata.
#'
#' 2. **Splitting ("split" method):**
#'    This method partitions the cells in each sample group into roughly equal,
#'    **non-overlapping subsets** (pseudo-replicates) by randomly splitting cells.
#'    It does not change count data or embeddings but simply annotates cells in metadata.
#'    Statistically, splitting provides multiple independent subsets for analysis,
#'    allowing evaluation of consistency or heterogeneity across subsets.
#'    This method is computationally lighter and preserves original embeddings (e.g., UMAP),
#'    since the counts are not altered.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param sample_column Character string naming the metadata column defining groups (e.g., "donor").
#' @param target_samples Character vector of group names to generate pseudo-replicates for. Defaults to all groups.
#' @param n_replicates Integer number of pseudo-replicates to generate per group.
#' @param method Character specifying the approach to create pseudo-replicates. Options:
#'   - `"bootstrap"`: Generate bootstrap samples with replacement, returns new Seurat object.
#'   - `"split"`: Partition cells into non-overlapping groups, updates metadata only.
#' @param seed Integer random seed for reproducibility.
#'
#' @return A Seurat object:
#' - If `method = "split"`, returns the input object with updated metadata column `pseudo_replicate`.
#' - If `method = "bootstrap"`, returns a new Seurat object containing bootstrap-sampled counts and metadata.
#'
#' @examples
#' # Split method - preserves embeddings, annotates metadata - Assess consistency across independent subsets.
#' obj_split <- generate_pseudo_replicates(obj, sample_column = "sample", n_replicates = 3, method = "split")
#'
#' # Bootstrap method - new object with resampled counts - Estimate variability & robustness of DE analysis and Perform rigorous statistical inference.
#' obj_boot <- generate_pseudo_replicates(obj, sample_column = "sample", n_replicates = 3, method = "bootstrap")
#'
#' Additional notes:
#'	•	Bootstrap is generally more statistically sound for quantifying uncertainty in DE but more computationally expensive.
#'  •	Split is often used as a quick diagnostic to ensure signals are not driven by a subset of cells.
#'  •	In practice, many analyses combine both: use split to identify consistency, then bootstrap for inference.
#' @export
generate_pseudo_replicates <- function(seurat_obj,
                                       sample_column = "donor",
                                       target_samples = NULL,
                                       n_replicates = 3,
                                       method = c("split", "bootstrap"),
                                       seed = 123) {
  library(Seurat)
  set.seed(seed)
  
  method <- match.arg(method)
  
  meta <- seurat_obj@meta.data
  if (!(sample_column %in% colnames(meta))) stop(paste("Column", sample_column, "not found in metadata"))
  meta[[sample_column]] <- as.character(meta[[sample_column]])
  
  if (is.null(target_samples)) target_samples <- unique(meta[[sample_column]])
  
  if (method == "split") {
    meta$pseudo_replicate <- meta[[sample_column]]  # initialize
    
    for (sample in target_samples) {
      sel <- which(meta[[sample_column]] == sample)
      n_cells <- length(sel)
      
      if (n_cells < n_replicates) {
        warning(sprintf("Sample %s has fewer cells (%d) than requested replicates (%d). Adjusting to %d replicates.",
                        sample, n_cells, n_replicates, n_cells))
        actual_reps <- n_cells
      } else {
        actual_reps <- n_replicates
      }
      
      shuffled_cells <- sample(sel)
      split_indices <- split(shuffled_cells, rep(1:actual_reps, length.out = length(shuffled_cells)))
      
      for (rep in seq_len(actual_reps)) {
        meta$pseudo_replicate[split_indices[[rep]]] <- paste0(sample, "_rep_", rep)
      }
    }
    
    seurat_obj@meta.data <- meta
    return(seurat_obj)
  }
  
  if (method == "bootstrap") {
    counts_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    all_cells_in_counts <- colnames(counts_data)
    
    all_cells <- list()
    all_meta <- list()
    
    for (sample in target_samples) {
      sel <- rownames(meta)[meta[[sample_column]] == sample]
      sel <- intersect(sel, all_cells_in_counts)
      n_cells <- length(sel)
      if (n_cells == 0) next
      
      for (rep in seq_len(n_replicates)) {
        boot_idx <- sample(sel, size = n_cells, replace = TRUE)
        rep_label <- paste0(sample, "_boot_", rep)
        
        counts_mat <- counts_data[, boot_idx, drop = FALSE]
        colnames(counts_mat) <- paste0(rep_label, "_", seq_along(boot_idx))
        
        meta_rep <- meta[boot_idx, , drop = FALSE]
        rownames(meta_rep) <- colnames(counts_mat)
        meta_rep$pseudo_replicate <- rep_label
        
        all_cells[[rep_label]] <- counts_mat
        all_meta[[rep_label]] <- meta_rep
      }
    }
    
    counts_combined <- do.call(cbind, all_cells)
    meta_combined <- do.call(rbind, all_meta)
    
    new_obj <- CreateSeuratObject(counts = counts_combined, meta.data = meta_combined)
    return(new_obj)
  }
}


EnrichGO <- function(de_table, 
                     gene_col = "gene",
                     cluster_col = "cell_ype",  
                     logFC_col = "avg_logFC",  
                     pval_col = "p_val_adj",  
                     pval_cutoff = 0.05, 
                     p_adj_method = "none",
                     ontology = "BP",  
                     min_genes = 10,  
                     output_dir = "Enrich_GO") {
  
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(dplyr)
  
  # Ensure required columns exist
  required_cols <- c(gene_col, cluster_col, logFC_col, pval_col)
  if (!all(required_cols %in% colnames(de_table))) {
    stop("Error: DE table must contain the specified columns: ", paste(required_cols, collapse=", "))
  }
  
  # Filter for significant genes
  de_table <- de_table %>% filter(.data[[pval_col]] < pval_cutoff)
  
  # Get unique clusters (e.g., cell types)
  unique_clusters <- unique(de_table[[cluster_col]])
  
  # Store results
  dir <- file.path(output_dir)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  go_results_list <- list()
  
  for (cluster in unique_clusters) {
    cat("Running GO enrichment for:", cluster,"...\n")
    
    # Separate Up and Down genes
    up_genes <- de_table %>%
      filter(.data[[cluster_col]] == cluster, .data[[logFC_col]] > 0) %>%
      pull(.data[[gene_col]])
    
    down_genes <- de_table %>%
      filter(.data[[cluster_col]] == cluster, .data[[logFC_col]] < 0) %>%
      pull(.data[[gene_col]])
    
    # Function to perform GO enrichment
    run_go <- function(gene_list, regulation_status) {
      if (length(gene_list) < min_genes) return(NULL)  # Skip if too few genes
      
      entrez_genes <- mapIds(org.Hs.eg.db, keys = gene_list, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
      entrez_genes <- na.omit(entrez_genes)  # Remove NAs
      
      if (length(entrez_genes) >= min_genes) {
        gse <- enrichGO(gene = entrez_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = ontology,  
                        pAdjustMethod = p_adj_method,
                        pvalueCutoff = pval_cutoff,
                        #qvalueCutoff = pval_cutoff,
                        readable = TRUE)
        
        gse@result$Description <- stringr::str_to_upper(gse@result$Description)
        go_result <- gse@result
        
        if (!is.null(go_result) && nrow(go_result) > 0) {
          go_df <- as.data.frame(go_result)
          
          # Only add metadata if results exist
          if (nrow(go_df) > 0) {
            go_df[[cluster_col]] <- cluster  # Add cluster column
            go_df$Regulation <- regulation_status  # "Up" or "Down"
            saveRDS(gse, file.path(dir, paste0(cluster, "_Enrich_GO_", ontology, ".rds")))
          }
          
          return(go_df)
        }
      }
      return(NULL)
    }
    
    # Run GO for Up and Down genes
    up_results <- run_go(up_genes, "Up")
    down_results <- run_go(down_genes, "Down")
    
    # Store only non-empty results
    if (!is.null(up_results) && nrow(up_results) > 0) {
      go_results_list[[paste0(cluster, "_Up")]] <- up_results
    }
    if (!is.null(down_results) && nrow(down_results) > 0) {
      go_results_list[[paste0(cluster, "_Down")]] <- down_results
    }
  }
  
  
  # Merge all results
  if (length(go_results_list) > 0) {
    final_go_results <- bind_rows(go_results_list)
    write.table(final_go_results, file.path(dir, paste0("Enrich_GO_", ontology, ".tsv")),sep = "\t", row.names = F, quote = F)
    cat("GO enrichment completed. Results saved in: ", dir,"\n")
    return(final_go_results)
  } else {
    cat("No significant GO terms found for any cell type.\n")
    return(NULL)
  }
}


#' Plot permutation test results for single-cell analysis
#'
#' This function creates a plot of observed log2 fold differences (log2FD) 
#' from a permutation test, highlighting statistically significant results.
#'
#' @param sc_utils_obj An object containing single-cell analysis results, 
#'   with a `@results$permutation` slot storing a `data.table` that includes:
#'   `FDR`, `obs_log2FD`, `boot_CI_2.5`, `boot_CI_97.5`, and `clusters`.
#' @param FDR_threshold Numeric, false discovery rate threshold for significance. Default is `0.05`.
#' @param log2FD_threshold Numeric, log2 fold-difference threshold for significance. 
#'   Default is `log2(1.5)`.
#' @param order_clusters Logical, if `TRUE` (default) orders clusters by observed log2FD.
#' @param theme.type Character string indicating the ggplot2 theme type to apply 
#'   via `plot_theme()`. Default is `"minimal"`.
#' @param font.size Numeric, base font size for the plot theme. Default is `10`.
#' @param legend.position Position of the legend (`"right"`, `"bottom"`, `"top"`, `"left"`). 
#'   Default is `"right"`.
#' @param legend.direction Direction of the legend (`"vertical"` or `"horizontal"`). 
#'   Default is `"vertical"`.
#' @param ... Additional arguments passed to `plot_theme()`.
#'
#' @return A `ggplot` object showing permutation results with confidence intervals 
#'   and significance highlighting.
#'
#' @details
#' The function:
#' 1. Copies the permutation results from the `sc_utils_obj` object.
#' 2. Marks results as significant if they pass both the FDR and log2FD thresholds.
#' 3. Optionally orders clusters by observed log2FD.
#' 4. Creates a point-range plot with bootstrapped confidence intervals and 
#'    color-coded significance.
#'
#' @examples
#' \dontrun{
#' p <- permut_plot(my_sc_utils_object, FDR_threshold = 0.01, log2FD_threshold = log2(2))
#' print(p)
#' }
#'
#' @export
permut_plot <- function(
    sc_utils_obj,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    order_clusters = TRUE,
    theme = "minimal",
    colors = c("salmon", "grey60"),
    font.size = 10,
    leg.ttl = 10,
    leg.size = 10,
    leg.pos = "right",
    leg.dir = "vertical",
    ...
) {
  ## Retrieve results.
  plot_data <- data.table::copy(sc_utils_obj@results$permutation)
  
  ## Mark the significant results.
  plot_data[, significance := ifelse(
    FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  )]
  
  plot_data[, significance := factor(significance, levels = c(
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  ))]
  
  ## Order the clusters by observed log2FD if requested.
  if (order_clusters) {
    plot_data[, clusters := forcats::fct_reorder(factor(clusters), dplyr::desc(obs_log2FD))]
  }
  
  ## Plot the results.
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
    plot_theme(
      theme.style = theme,
      font.size = font.size,
      leg.pos = leg.pos,
      leg.dir = leg.dir,
      leg.size = leg.size,
      y.ttl = F,
      ...
    ) +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = colors) +
    coord_flip()
  
  return(p)
}


#' Find Marker Genes in a Seurat Object
#' This function identifies marker genes for specified clusters in a Seurat object
#' using differential expression testing within each cluster.
#' Run Seurat::FindMarkers across clusters with multiple tests + consensus support + visualization
#'
#' @param object Seurat object
#' @param group.by Metadata column for clustering/identity
#' @param test.use Differential test(s) ("wilcox", "bimod", "roc", "t",
#'   "negbinom", "poisson", "LR", "MAST" or "all")
#' @param only.pos Keep only positive markers (default TRUE)
#' @param min.pct Minimum fraction of cells expressing a gene (default 0.25)
#' @param min.diff.pct Minimum difference in pct (default -Inf)
#' @param man.logfc.threshold LogFC threshold (default 0.25)
#' @param clusters.to.exclude Clusters to skip (default none)
#' @param max.cells.per.ident Downsample max cells per cluster (default Inf)
#' @param consensus Logical, if TRUE build consensus markers across tests
#' @param consensus_min_tests Integer, minimum number of tests a gene must be significant in
#' @param alpha Adjusted p-value threshold for significance (default 0.05)
#' @param return_both Logical, return both raw + consensus results (default FALSE)
#' @param plot_type "none", "bar", or "upset" (default "none")
#' @param plot_cluster Cluster name for upset plot (default first cluster)
#' @param ... Passed to Seurat::FindMarkers
#'
#' @return Data frame(s) of markers, optionally with ggplot object
#' @export
cellmarker <- function(object,
                       group.by,
                       assay = "RNA",
                       features = NULL,
                       test.use = 'wilcox',
                       only.pos = TRUE,
                       min.pct = 0.01,
                       min.diff.pct = -Inf,
                       logfc.threshold = 0.1,
                       clusters.to.exclude = c(),
                       max.cells.per.ident = Inf, 
                       latent.vars = NULL,
                       consensus = FALSE,
                       consensus_min_tests = 2,
                       alpha = 0.05,
                       return_both = FALSE,
                       plot_type = c("none","bar","upset"),
                       plot_cluster = NULL,
                       ...) {
  
  requireNamespace("ggplot2")
  if ("upset" %in% plot_type) requireNamespace("UpSetR")
  
  if (!inherits(object, "Seurat")) stop("The provided object is not a Seurat object.")
  
  if (!missing(group.by)) object <- Seurat::SetIdent(object = object, value = group.by)
  
  valid_tests <- c("wilcox", "bimod", "roc", "t", 
                   "negbinom", "poisson", "LR", "MAST")
  
  if (identical(test.use, "all")) {
    tests_to_run <- valid_tests
  } else {
    if (!all(test.use %in% valid_tests)) stop(paste("Invalid test.use. Choose from:", paste(valid_tests, collapse = ", ")))
    tests_to_run <- test.use
  }
  
  clusters.to.test <- sort(unique(object@active.ident))
  clusters.to.test <- setdiff(clusters.to.test, clusters.to.exclude)
  
  results_all <- list()
  
  for (test in tests_to_run) {
    message("=== Running test: ", test, " ===")
    joined <- data.frame()
    
    for (i in seq_along(clusters.to.test)) {
      cluster <- clusters.to.test[i]
      message(sprintf("[ %d / %d ] %s - cluster: %s ...", 
                      i, length(clusters.to.test), test, cluster))
      
      markers <- Seurat::FindMarkers(
        object,
        ident.1 = cluster, 
        ident.2 = NULL, 
        features = features,
        only.pos = only.pos, 
        assay = assay, 
        slot = "data",
        test.use = test,
        min.pct = min.pct,
        min.cells.group = 3,
        min.diff.pct = min.diff.pct,
        logfc.threshold = logfc.threshold,
        max.cells.per.ident = max.cells.per.ident, 
        latent.vars = latent.vars,
        ...
      )
      
      if (!is.null(markers) && nrow(markers) > 0) {
        markers$cluster <- cluster
        markers$gene <- rownames(markers)
        markers$test.use <- test
        rownames(markers) <- NULL
        joined <- rbind(joined, markers)
      } else {
        message(sprintf("Skipping cluster '%s': no markers found.", cluster))
      }
    }
    
    results_all[[test]] <- joined
  }
  
  all_results <- dplyr::bind_rows(results_all)
  
  #- Single test shortcut-
  if (length(tests_to_run) == 1 && !consensus && !return_both) {
    return(all_results)
  }
  
  #- Consensus / visualization-
  consensus_df <- NULL
  if (consensus || return_both) {
    consensus_df <- all_results %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::group_by(cluster, gene) %>%
      dplyr::summarise(
        n_tests = dplyr::n_distinct(test.use),
        mean_logFC = mean(avg_log2FC, na.rm = TRUE),
        min_p_val_adj = min(p_val_adj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_tests >= consensus_min_tests)
  }
  
  plot_type <- match.arg(plot_type)
  plt <- NULL
  
  if (plot_type == "bar") {
    marker_counts <- all_results %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::group_by(cluster, test.use) %>%
      dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")
    
    plt <- ggplot2::ggplot(marker_counts, 
                           ggplot2::aes(x = cluster, y = n_genes, fill = test.use)) +
      ggplot2::geom_bar(stat="identity", position="dodge") +
      plot_theme(theme.type = "classic", x.angle = 45) +
      ggplot2::labs(title="Significant markers per test per cluster", 
                    y="Number of markers", x="Cluster")
  }
  
  if (plot_type == "upset") {
    if (is.null(plot_cluster)) plot_cluster <- clusters.to.test[1]
    df_upset <- all_results %>%
      dplyr::filter(cluster == plot_cluster, p_val_adj <= alpha) %>%
      dplyr::select(gene, test.use) %>%
      dplyr::distinct()
    
    mat <- table(df_upset$gene, df_upset$test.use) > 0
    plt <- UpSetR::upset(UpSetR::fromMatrix(mat), 
                         mainbar.y.label = paste("Overlap of markers -", plot_cluster))
  }
  
  if (return_both) return(list(raw = all_results, consensus = consensus_df, plot = plt))
  if (consensus) return(list(consensus = consensus_df, plot = plt))
  
  return(list(raw = all_results, plot = plt))
}


plot_de_summary <- function(
    DE, 
    cell.type.col = "celltype",
    pvalue.col = "p_val_adj",
    logFC.col = "avg_log2FC",
    pvalue = 0.05, 
    logFC = 0, 
    plot.title = NULL,
    font.size = 12,
    axis.text.size = 12,
    theme = "test",
    x.angle = 90,
    leg.pos = "right",
    leg.dir = "vertical",
    custom.colors = NULL,
    cell.type.order = NULL,
    flip.coords = FALSE,
    dodge.plot = FALSE,
    show.text = TRUE
) {
  
  # Harmonize column names: remove '_' and '.'
  colnames(DE) <- gsub("[._]", "", colnames(DE))
  cell.type.col <- gsub("[._]", "", cell.type.col)
  pvalue.col    <- gsub("[._]", "", pvalue.col)
  logFC.col     <- gsub("[._]", "", logFC.col)
  
  # REMOVE NA CELL TYPES HERE
  DE <- DE[!is.na(DE[[cell.type.col]]), ]
  
  # Check existence
  required <- c(cell.type.col, pvalue.col, logFC.col)
  if (!all(required %in% colnames(DE))) {
    stop("Missing column(s): ", paste(required[!required %in% colnames(DE)], collapse = ", "))
  }
  
  # Ensure no list columns (convert to character)
  if (is.list(DE[[cell.type.col]])) {
    DE[[cell.type.col]] <- as.character(DE[[cell.type.col]])
  }
  
  # Assign DE direction
  DE$updown <- dplyr::case_when(
    DE[[pvalue.col]] <= pvalue & DE[[logFC.col]] >  logFC ~ "Up",
    DE[[pvalue.col]] <= pvalue & DE[[logFC.col]] < -logFC ~ "Down",
    TRUE ~ "NS"
  )
  
  df.sig <- DE[DE$updown != "NS", ]
  
  # SAFE COUNTING (NO tidy-eval) 
  df.m <- df.sig %>%
    dplyr::group_by(.data[[cell.type.col]], updown) %>%
    dplyr::summarise(value = dplyr::n(), .groups = "drop")
  
  # Rename safely (direct, no tidy-eval)
  names(df.m)[names(df.m) == cell.type.col] <- "cell.type"
  
  # Filter positive
  df.m <- df.m[df.m$value > 0, ]
  
  # Order
  if (!is.null(cell.type.order)) {
    df.m$cell.type <- factor(df.m$cell.type, levels = cell.type.order)
  } else {
    df.m$cell.type <- factor(df.m$cell.type)
  }
  
  df.m$updown <- factor(df.m$updown, c("Up", "Down"))
  
  # Default colors
  if (is.null(custom.colors)) {
    custom.colors <- c("Up"="#B2182B","Down"="#2166AC")
  }
  
  # Plot
  if (dodge.plot) {
    
    p <- ggplot(df.m, aes(x = cell.type, y = value, fill = updown)) +
      geom_bar(stat="identity", position=position_dodge(width=0.8),
               color="black", linewidth=0.2) +
      scale_fill_manual(values = custom.colors)
    
    if (show.text) {
      p <- p +
        geom_text(aes(label=value, group=updown),
                  position=position_dodge(width=0.8),
                  vjust=-0.4, size=3)
    }
    
  } else {
    
    p <- ggplot(df.m, aes(x = cell.type)) +
      geom_bar(data=df.m[df.m$updown=="Up",],
               aes(y=value, fill="Up"),
               stat="identity", color="black",  width = 0.5, linewidth=0.2) +
      geom_bar(data=df.m[df.m$updown=="Down",],
               aes(y=-value, fill="Down"),
               stat="identity", color="black", width = 0.5, linewidth=0.2, alpha=0.8) +
      geom_hline(yintercept=0, color="grey80") +
      scale_fill_manual(values = custom.colors)
    
  }
  
  # Theme
  p <- p +
    xlab("") + ylab("Number of DE genes") +
    ggtitle(plot.title) +
    labs(fill="DE genes") +
    plot_theme(theme=theme, x.angle=x.angle,
               font.size=font.size, 
               leg.pos=leg.pos,
               leg.dir=leg.dir)
  
  if (flip.coords) p <- p + coord_flip()
  
  return(p)
}


#' Plot Differential Expression Results for Pseudobulk Data
#'
#' Generates volcano plots, heatmaps, or both for pseudobulk differential expression (DE) results,
#' highlighting significant genes based on user-defined thresholds.
#'
#' @param DE A data frame containing differential expression results. Must include at least:
#'   \code{cell_type}, \code{gene}, \code{avg_logFC}, and the column specified by \code{p_value_column}.
#' @param matrices A named list of count matrices (one per cell type). Row names are genes, column names are samples.
#' @param n_lab Integer. Maximum number of genes to label (split evenly between up- and down-regulated). Default is 100.
#' @param width Numeric. Width of output plots in inches. Default is 12.
#' @param height Numeric. Height of output plots in inches. Default is 6.
#' @param base_size Numeric. Base font size for plot text. Default is 12.
#' @param label_size Numeric. Font size for gene labels in volcano plots. Default is 8.
#' @param anno_legend Logical. Whether to display annotation legends in heatmaps. Default is \code{FALSE}.
#' @param legend Logical. Whether to display the main heatmap legend. Default is \code{TRUE}.
#' @param show.rownames Logical. Whether to display row (gene) names in heatmaps. Default is \code{FALSE}.
#' @param show.colnames Logical. Whether to display column (sample) names in heatmaps. Default is \code{FALSE}.
#' @param angle.col Numeric. Rotation angle for column labels in heatmaps. Default is 90.
#' @param pval_cutoff Numeric. Adjusted p-value threshold for significance. Default is 0.05.
#' @param logfc_cutoff Numeric. Absolute log fold change threshold for significance. Default is 0.
#' @param label_rectangle Logical. Currently unused. Default is \code{FALSE}.
#' @param min_samples_gene Integer. Minimum number of samples in which a gene must be expressed (|value| >= 1) in either group. Default is 3.
#' @param group1_pattern Character. Pattern to match column names for group 1 samples. Default is \code{"condition1"}.
#' @param group2_pattern Character. Pattern to match column names for group 2 samples. Default is \code{"condition2"}.
#' @param out_dir Character. Output directory where plots and gene lists will be saved.
#' @param border_color Color for heatmap cell borders. Default is \code{NA}.
#' @param dynamic_height Numeric or \code{NULL}. Custom heatmap height in inches. If \code{NULL}, calculated dynamically based on number of genes. Default is \code{NULL}.
#' @param highlight_genes Optional character vector of gene names to highlight in heatmaps.
#' @param custom_annotation_colors Optional list specifying colors for heatmap annotations.
#' @param color_theme Character. Name of viridis color palette option (see \code{\link[viridis]{viridis}}). Default is \code{"viridis"}.
#' @param plot_type Character. Type of plot to generate: \code{"heatmap"}, \code{"volcano"}, or \code{"both"}. Default is \code{"both"}.
#' @param p_value_column Character. Column name in \code{DE} containing the p-values or adjusted p-values to filter by. Default is \code{"p_val_adj"}.
#'
#' @details
#' For each \code{cell_type} in \code{DE}:
#' \enumerate{
#'   \item Filters genes based on expression in each group and the specified significance thresholds.
#'   \item Selects the top \code{n_lab} significant genes (half up-, half down-regulated).
#'   \item Generates:
#'     \itemize{
#'       \item A volcano plot with optional gene labels.
#'       \item A heatmap of significant genes (rows) by samples (columns).
#'     }
#'   \item Saves plots to \code{out_dir} in PDF format and the list of significant genes as an RDS file.
#' }
#'
#' Volcano plots use \code{ggplot2} and \code{ggrepel} for labeled points.
#' Heatmaps use \code{pheatmap} with row scaling and optional annotation colors.
#'
#' @return
#' Saves plots and gene lists to the specified output directory. Invisibly returns \code{NULL}.
#'
#' @import ggplot2 dplyr cowplot viridis ggrepel pheatmap edgeR
#' @importFrom stringr str_split_fixed
#' @export
#'
#' @examples
#' \dontrun{
#' Plot_DE_Pseudobulk(
#'   DE = de_results,
#'   matrices = count_matrices,
#'   out_dir = "plots/",
#'   pval_cutoff = 0.01,
#'   logfc_cutoff = 1,
#'   plot_type = "both"
#' )
#' }
#' @export
#' 
Plot_DE_Pseudobulk <- function(DE, 
                               matrices, 
                               n_lab = 100, 
                               width = 12, 
                               height = 6, 
                               base_size = 12, 
                               label_size = 8,
                               anno_legend = F, 
                               legend = T, 
                               show.rownames = F, 
                               show.colnames = F, 
                               angle.col = 90,
                               pval_cutoff = 0.05, 
                               logfc_cutoff = 0, 
                               label_rectangle = F, 
                               min_samples_gene = 3,
                               group1_pattern = "condition1", 
                               group2_pattern = "condition2", 
                               out_dir, 
                               border_color = NA,
                               dynamic_height = NULL,
                               highlight_genes = NULL,
                               custom_annotation_colors = NULL, 
                               color_theme = "viridis",
                               plot_type = "both", 
                               p_value_column = "p_val_adj") {
  
  cell_types <- unique(DE$cell_type)
  
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    DE_sub <- DE[DE$cell_type == cell_type, ]
    
    cat("Processing cell type:", cell_type, "\n")
    cat("Total genes in DE_sub:", nrow(DE_sub), "\n")
    
    # Select genes for groups
    mtx_group2 <- dplyr::select(matrices[[cell_type]], contains(group2_pattern))
    mtx_group2 <- mtx_group2[apply(abs(mtx_group2) >= 1, 1, function(x) sum(x) >= min_samples_gene), ]
    group2_genes <- rownames(mtx_group2)
    
    mtx_group1 <- dplyr::select(matrices[[cell_type]], contains(group1_pattern))
    mtx_group1 <- mtx_group1[apply(abs(mtx_group1) >= 1, 1, function(x) sum(x) >= min_samples_gene), ]
    group1_genes <- rownames(mtx_group1)
    
    # Filter significant genes
    sig <- dplyr::filter(DE_sub, .data[[p_value_column]] <= pval_cutoff & gene %in% c(group1_genes, group2_genes))
    
    # Split DE_sub into up/down regulated
    DE_pos <- dplyr::filter(DE_sub, .data[[p_value_column]] <= pval_cutoff & avg_logFC > logfc_cutoff) %>% 
      dplyr::arrange(desc(avg_logFC)) %>%
      dplyr::top_n(n = ceiling(n_lab / 2), wt = avg_logFC) %>% dplyr::pull(gene)
    DE_pos <- intersect(DE_pos, sig$gene)
    
    DE_neg <- dplyr::filter(DE_sub, .data[[p_value_column]] <= pval_cutoff & avg_logFC < logfc_cutoff) %>% 
      dplyr::arrange(avg_logFC) %>% 
      dplyr::top_n(n = floor(n_lab / 2), wt = avg_logFC) %>% dplyr::pull(gene)
    DE_neg <- intersect(DE_neg, sig$gene)
    
    # Log-normalized counts
    norm_counts <- log(edgeR::cpm(matrices[[cell_type]]) + 1, 2)
    colnames(norm_counts) <- sub(":*", "", colnames(norm_counts))
    
    # Highlight genes: include only if significant in DE_sub
    sig_highlight <- NULL
    if (!is.null(highlight_genes)) {
      sig_highlight <- highlight_genes[highlight_genes %in% DE_sub$gene &
                                         DE_sub[[p_value_column]][match(highlight_genes, DE_sub$gene)] <= pval_cutoff]
      if (length(sig_highlight) > 0) {
        message("Significant highlight genes for ", cell_type, ": ", paste(sig_highlight, collapse = ", "))
      }
    }
    
    # Combine DE_pos, DE_neg, and significant highlight genes
    sig_genes <- unique(c(DE_pos, DE_neg, sig_highlight))
    
    # Ensure genes exist in the matrix
    sig_genes <- sig_genes[sig_genes %in% rownames(norm_counts)]
    
    # Warn if significant highlight genes are missing from matrix
    missing_genes <- setdiff(sig_highlight, rownames(norm_counts))
    if (length(missing_genes) > 0) {
      warning("The following significant highlight genes are not in the matrix and will be skipped: ", 
              paste(missing_genes, collapse = ", "))
    }
    
    if (length(sig_genes) > 1) {
      
      safe_cell_type <- gsub("_", " ", cell_type)
      
      # Volcano plot
      df <- rowMeans(norm_counts)
      df <- data.frame(baseMean = df, gene = rownames(norm_counts))
      merged <- merge(DE_sub, df, by = "gene")
      volcanoplot <- merged[, c("gene", "baseMean", "avg_logFC", "p_val_adj")]
      colnames(volcanoplot) <- c("gene", "baseMean", "log2FoldChange", "padj")
      
      volcanoplot$Include <- with(volcanoplot, ifelse(padj <= pval_cutoff & !(gene %in% sig_genes), FALSE, TRUE))
      volcanoplot$log_padj <- -log10(volcanoplot$padj)
      volcanoplot$Significance <- with(volcanoplot, ifelse(padj <= pval_cutoff & abs(log2FoldChange) >= logfc_cutoff,
                                                           ifelse(log2FoldChange > 0, "Up", "Down"), "NS"))
      
      label_selected <- subset(volcanoplot, Include & Significance != "NS")
      
      volcano <- ggplot(subset(volcanoplot, Include), aes(x = log2FoldChange, y = log_padj)) +
        geom_point(aes(color = Significance), stroke = 0.5, size = 1.5) +
        scale_color_manual(values = c("Up" = "#DD5D74", "Down" = "#269BCF", "NS" = "gray")) +
        theme_bw() +
        theme(plot.title = element_text(size = base_size, face = "bold", hjust = 0.5),
              axis.title = element_text(size = base_size),
              axis.text = element_text(size = base_size),
              legend.text = element_text(size = base_size),
              panel.grid = element_blank(),
              legend.position = "right",
              legend.title = element_blank()) +
        labs(x = bquote(~Log[2]~"FC ["*.(group1_pattern)*"/"*.(group2_pattern)*"]"), 
             y = bquote(~-Log[10]~italic("FDR"))) +
        geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
        geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed")
      
      if (nrow(label_selected) > 0) {
        volcano <- volcano + ggrepel::geom_text_repel(
          data = label_selected,
          aes(label = gene, fontface = "italic"),
          size = label_size, color = "black", bg.color = "white",
          bg.r = 0.2, box.padding = 0.35, point.padding = 0.5,
          segment.color = 'grey50'
        )
      }
      
      # Heatmap
      meta <- as.data.frame(str_split_fixed(colnames(norm_counts), ":", 2))
      rownames(meta) <- colnames(norm_counts)
      colnames(meta) <- c("Samples", "Groups")
      mat <- norm_counts[unique(sig_genes), ]
      mat <- mat[rowSums(mat != 0) > 0, ]
      mat <- mat[apply(mat, 1, var) > 0, ]
      
      newnames <- sapply(rownames(mat), function(g) {
        if (!is.null(sig_highlight) && g %in% sig_highlight) {
          as.expression(bquote(italic(.(g)) ~ bold("*")))
        } else {
          as.expression(bquote(italic(.(g))))
        }
      })
      
      if (is.null(dynamic_height)) {
        num_genes <- nrow(mat)
        dynamic_height <- max(6, min(0.2 * num_genes, 20))
      }
      
      heatmap <- pheatmap::pheatmap(
        mat = mat, annotation_col = meta, annotation_colors = custom_annotation_colors, 
        clustering_method = "complete", treeheight_row = 0, cluster_row = T, 
        cluster_col = F, angle_col = angle.col, show_rownames = show.rownames, 
        show_colnames = show.colnames, scale = 'row', legend = legend, 
        annotation_legend = anno_legend, fontsize = base_size, fontsize_col = base_size,
        fontsize_row = base_size, legend_breaks = seq(-2, 2, by = 1), legend_labels = seq(-2, 2, by = 1),
        border_color = border_color, silent = T, main = paste(cell_type),
        color = viridis::viridis(100, option = color_theme),
        breaks = seq(-2, 2, length.out = 101), labels_row = as.expression(newnames)
      )
      
      # Save significant genes
      output_path <- out_dir
      if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
      saveRDS(sig_genes, file.path(output_path, paste0(cell_type, "_sig_genes.rds")))
      
      # Save plots
      if (plot_type %in% c("heatmap", "both")) {
        tiff(file.path(output_path, paste0(cell_type, "__heatmap.tiff")),
             width = width, height = dynamic_height, units = "in", res = 300)
        grid::grid.newpage()
        grid::grid.draw(heatmap[[4]])
        dev.off()
        
        pdf(file.path(output_path, paste0(cell_type, "__heatmap.pdf")), width = width, height = dynamic_height)
        grid::grid.newpage()
        grid::grid.draw(heatmap[[4]])
        dev.off()
      }
      
      if (plot_type %in% c("volcano", "both")) {
        pdf(file.path(output_path, paste0(cell_type, "_volcano.pdf")), width = width, height = height)
        print(volcano + ggtitle(safe_cell_type))
        dev.off()
      }
      
      if (plot_type == "both") {
        combined_plot <- cowplot::plot_grid(volcano, heatmap[[4]], ncol = 2, rel_widths = c(1, 1))
        pdf(file.path(output_path, paste0(cell_type, "_volcano_heatmap.pdf")), width = width, height = dynamic_height)
        print(combined_plot)
        dev.off()
      }
      
    }
  }
}


#' Pseudobulk differential expression using Libra with explicit reference control
#'
#' Performs pseudobulk differential expression analysis across cell types using
#' \pkg{Libra}, with explicit control of the reference and test conditions.
#' Log fold-changes are oriented such that positive values indicate upregulation
#' in the test condition relative to the reference.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Aggregates counts per cell type and sample (pseudobulk)
#'   \item Filters cell types with insufficient cells per condition
#'   \item Runs differential expression using \code{Libra::run_de()}
#'   \item Applies expression and sample-level gene filtering
#'   \item Extracts significant genes and balanced top-N markers per cell type
#'   \item Returns DE tables, pseudobulk matrices, summaries, and joined outputs
#' }
#'
#' @param object A \code{\link[Seurat]{Seurat}} object.
#' @param ident_col Metadata column used to define cell identities.
#' @param sample_col Metadata column defining biological replicates.
#' @param label_col Metadata column defining experimental condition.
#' @param ref_level Character scalar specifying the reference (baseline) condition.
#' @param test_level Character scalar specifying the test condition.
#' @param assay Assay name used for DE analysis. Default is \code{"RNA"}.
#' @param slot Slot name containing count data. Default is \code{"counts"}.
#' @param de_method Differential expression method passed to \pkg{Libra}.
#'   Default is \code{"edgeR"}.
#' @param de_type Type of test for DE analysis (e.g. \code{"LRT"}, \code{"QLF"}).
#' @param min_cells Minimum number of cells per condition required for a cell type
#'   to be tested.
#' @param min_samples_gene Minimum number of pseudobulk samples in which a gene
#'   must be expressed to be considered.
#' @param min_expression Minimum expression threshold used for gene filtering.
#' @param logfc_cutoff Absolute log fold-change cutoff for significance.
#' @param fdr_cutoff Adjusted p-value (FDR) cutoff for significance.
#' @param n_lab Number of genes per cell type retained for plotting (balanced
#'   between up- and down-regulated).
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{contrast}{Character string describing the contrast (test vs reference).}
#'   \item{DE}{Full differential expression table from \pkg{Libra}.}
#'   \item{DE_sig_all}{All significant genes passing thresholds.}
#'   \item{DE_sig_plot}{Balanced top-N genes per cell type for plotting.}
#'   \item{matrices}{Pseudobulk expression matrices per cell type.}
#'   \item{meta}{Metadata used for DE analysis.}
#'   \item{joined}{Pseudobulk matrices joined with DE statistics.}
#'   \item{upset_genes}{List of significant genes per cell type.}
#'   \item{DE_summary}{Summary table of up/down-regulated genes per cell type.}
#' }
#'
#' @note
#' Positive \code{avg_logFC} values indicate upregulation in the \code{test_level}
#' relative to the \code{ref_level}.
#'
#' @examples
#' \dontrun{
#' pb <- Run_DE_Pseudobulk(
#'   object = seurat_obj,
#'   ident_col = "celltype",
#'   sample_col = "sample",
#'   label_col = "group",
#'   ref_level = "CTRL",
#'   test_level = "CF"
#' )
#' }
#'
#' @import Seurat
#' @import Libra
#' @import dplyr
#' @import tibble
#' @import plyr
#'
#' @export
#' 
cellbulkde <- function(
    object,
    ident_col,
    sample_col,
    label_col,
    ref_level,
    test_level,
    assay = "RNA",
    slot = "counts",
    de_method = "edgeR",
    de_type = "LRT",
    min_cells_per_group = 10,
    min_samples_gene = 3,
    min_expression = 1,
    logfc_cutoff = 0.5,
    fdr_cutoff = 0.05
) {
  
  stopifnot(inherits(object, "Seurat"))
  
  ## Set identities 
  Idents(object) <- ident_col
  
  ## Build metadata 
  meta <- object@meta.data %>%
    transmute(
      cell_type = factor(Idents(object)),
      replicate = .data[[sample_col]],
      label = factor(.data[[label_col]], levels = c(ref_level, test_level))
    )
  
  ## Filter cell types by minimum cells
  valid_cell_types <- meta %>%
    dplyr::count(cell_type, label) %>%
    tidyr::pivot_wider(names_from = label, values_from = n, values_fill = 0) %>%
    dplyr::filter(if_all(all_of(c(ref_level, test_level)),
                         ~ . >= min_cells_per_group)) %>%
    pull(cell_type)
  
  if (length(valid_cell_types) == 0)
    stop("No cell types meet minimum cell requirements.")
  
  ## Subset object + metadata
  object <- subset(object, idents = valid_cell_types)
  meta <- meta[colnames(object), ]
  
  ## Counts
  counts <- GetAssayData(object, assay = assay, layer = slot)
  gc()
  
  ## Run DE
  deg <- Libra::run_de(
    input = counts,
    meta = meta,
    de_family = "pseudobulk",
    de_method = de_method,
    de_type = de_type
  )
  
  deg$contrast <- paste(test_level, "vs", ref_level)
  
  ## Pseudobulk matrices
  matrices <- Libra::to_pseudobulk(counts, meta)
  
  # Per cell type processing 
  deg_sig <- upset <- summary <- list()
  
  # Per cell type processing
  deg_sig <- res_plot <- upset <- summary <- list()
  
  for (ct in unique(deg$cell_type)) {
    
    deg_ct <- dplyr::filter(deg, cell_type == ct)
    
    # Expression filter
    keep_genes <- rownames(matrices[[ct]])[
      rowSums(abs(matrices[[ct]]) >= min_expression) >= min_samples_gene
    ]
    
    sig <- deg_ct %>%
      dplyr::filter(
        gene %in% keep_genes,
        p_val_adj <= fdr_cutoff,
        abs(avg_logFC) >= logfc_cutoff
      )
    
    deg_sig[[ct]] <- sig
    upset[[ct]]   <- sig$gene
    
  }
  
  # Join pseudobulk + deg
  joined <- plyr::rbind.fill(lapply(names(matrices), function(ct) {
    m <- matrices[[ct]][deg$gene[deg$cell_type == ct], , drop = FALSE]
    m <- tibble::rownames_to_column(as.data.frame(m), "gene")
    merge(m, filter(deg, cell_type == ct), by = "gene")
  }))
  joined[is.na(joined)] <- 0
  
  joined <- joined %>%
    dplyr::select(cell_type, gene, avg_logFC, p_val, p_val_adj, de_family, de_method, de_type, contrast, everything())
  
  return(list(
    contrast = paste(test_level, "vs", ref_level),
    deg = deg,
    joined = joined,
    deg_sig = bind_rows(deg_sig),
    matrices = matrices,
    upset_genes = upset
  ))
}


plot_go_terms <- function(go_df,
                          go_terms_of_interest = NULL,
                          cluster,
                          top_n = 10,
                          wrap_width = 50,
                          color_values = c("Up" = "#740001", "Down" = "#6497b1"),
                          size_range = c(2, 6),
                          theme = "bw",
                          x_angle = 0,
                          font_size = 12,
                          plot_title_size = 16) {
  
  go_df$Description <- stringr::str_to_upper(go_df$Description)
  
  # Filter for the target cell type
  df <- go_df %>%
    dplyr::filter(.data[["cell_type"]] == cluster)
  
  # If terms are NULL, select top_n by FoldEnrichment
  if (is.null(go_terms_of_interest)) {
    go_terms_of_interest <- df %>%
      dplyr::arrange(desc(.data[["FoldEnrichment"]])) %>%
      dplyr::pull(.data[["Description"]]) %>%
      unique() %>%
      head(top_n)
  }
  
  # Filter for selected terms
  df <- df %>%
    dplyr::filter(.data[["Description"]] %in% go_terms_of_interest) %>%
    dplyr::arrange(desc(.data[["FoldEnrichment"]])) %>%
    dplyr::group_by(.data[["Description"]]) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      log10FDR = -log10(.data[["pvalue"]]),
      Count = as.integer(.data[["Count"]]),
      CountSize = Count * 0.8,
      Regulation = factor(.data[["Regulation"]], levels = c("Up", "Down"))
    )
  
  # Wrap long GO term names
  df[["Description"]] <- stringr::str_wrap(df[["Description"]], width = wrap_width)
  
  # Plot
  p <- ggplot(df, aes(x = log10FDR, y = reorder(.data[["Description"]], log10FDR))) +
    geom_segment(aes(x = 0, xend = log10FDR,
                     yend = reorder(.data[["Description"]], log10FDR),
                     color = Regulation),
                 size = 0.8, alpha = 1) +
    geom_point(aes(size = .data[["FoldEnrichment"]], color = Regulation), alpha = 1) +
    scale_size_identity(name = "Gene Count", guide = "legend",
                        breaks = df[["FoldEnrichment"]], labels = df$Count) +
    scale_color_manual(values = color_values) +
    scale_size(range = size_range) +
    labs(
      x = expression(-log[10](FDR)),
      y = "",
      size = "Enrichment",
      color = "Regulation",
      title = cluster
    ) +
    plot_theme(theme = theme, x.angle = x_angle, font.size = font_size)
  
  return(p)
}


#' Subset a Seurat object with optional preprocessing and v5-safe layer repair
#'
#' This function performs a robust subset of a Seurat object, automatically
#' repairing corrupted assay layers commonly encountered in Seurat v5.
#' Optionally, it runs a full preprocessing pipeline including normalization,
#' variable feature selection, scaling, dimensional reduction (PCA or Harmony),
#' clustering, and UMAP embedding.
#'
#' The function unifies simple identity-based subsetting and metadata-based
#' filtering with a safe, retryable mechanism that prevents failures due to
#' malformed assay layers.
#'
#' @param object A \code{Seurat} object.
#' @param subset.expr An optional logical expression for subsetting cells based
#'   on metadata (passed to \code{\link[Seurat]{subset}} via \code{subset=}).
#' @param idents Optional identity class values to keep or remove.
#' @param invert Logical; whether to invert the identity selection.
#' @param group.by Optional metadata column used to set identities prior to
#'   subsetting.
#' @param preprocess Logical; whether to run preprocessing on the subset
#'   (default \code{TRUE}).
#' @param use.harmony Logical; whether to use Harmony for dimensional reduction
#'   (default \code{TRUE}).
#' @param harmony.group Metadata column used for Harmony batch correction
#'   (default \code{"orig.ident"}).
#' @param ndims Number of dimensions (PCs) to use for PCA/Harmony, neighbors,
#'   clustering, and UMAP (default \code{50}).
#' @param resolution Clustering resolution for
#'   \code{\link[Seurat]{FindClusters}} (default \code{0.5}).
#' @param min.dist Minimum distance parameter for UMAP (default \code{0.3}).
#' @param spread Spread parameter for UMAP (default \code{1}).
#'
#' @return A subsetted \code{Seurat} object. If \code{preprocess = TRUE}, the
#'   object is normalized, scaled (if applicable), reduced (PCA or Harmony),
#'   clustered, and embedded with UMAP.
#'
#' @details
#' This function provides:
#' \itemize{
#'   \item Safe subsetting compatible with Seurat v5 layer architecture.
#'   \item Automatic removal of corrupted or incompatible assay layers.
#'   \item Identity- or metadata-based subsetting.
#'   \item Optional preprocessing and dimensional reduction.
#'   \item Automatic fallback to PCA when Harmony is not applicable.
#' }
#'
#' Scaling is skipped when only a single identity level is present.
#'
#' @examples
#' \dontrun{
#' # Simple identity-based subset
#' sub <- subset_seurat2(obj, idents = "Mesenchyme")
#'
#' # Metadata-based subset without preprocessing
#' sub <- subset_seurat2(
#'   obj,
#'   subset.expr = ann1 == "Mesenchyme" & nCount_RNA > 1000,
#'   preprocess = FALSE
#' )
#'
#' # Full preprocessing with Harmony
#' sub <- subset_seurat2(
#'   obj,
#'   group.by = "ann2",
#'   idents = c("A", "B"),
#'   preprocess = TRUE,
#'   use.harmony = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[Seurat]{subset}},
#' \code{\link[Seurat]{NormalizeData}},
#' \code{\link[harmony]{RunHarmony}}
#'
#' @export
#' 
cellslice <- function(
    object,
    subset.expr = NULL,
    idents = NULL,
    invert = FALSE,
    group.by = NULL,
    preprocess = TRUE,
    use.harmony = TRUE,
    harmony.group = "orig.ident",
    ndims = 50,
    resolution = 0.5,
    run.umap = TRUE,
    min.dist = 0.3,
    spread = 1
) {
  stopifnot(inherits(object, "Seurat"))
  
  # Identity handling
  if (!is.null(group.by)) {
    object <- SetIdent(object, value = group.by)
  }
  
  DefaultAssay(object) <- "RNA"
  
  # Safe subset helper (Seurat v5 layer-safe)
  safe_subset <- function(obj) {
    tryCatch(
      {
        if (!is.null(idents)) {
          subset(obj, idents = idents, invert = invert)
        } else if (!is.null(subset.expr)) {
          subset(obj, subset = subset.expr)
        } else {
          obj
        }
      },
      error = function(e) {
        if (grepl("incorrect number of dimensions", e$message)) {
          return(NULL)
        }
        stop(e)
      }
    )
  }
  
  # First attempt
  result <- safe_subset(object)
  
  # Repair corrupted assay layers if needed
  if (is.null(result)) {
    warning("Subsetting failed. Repairing corrupted assay layers...", call. = FALSE)
    
    for (assay_name in names(object@assays)) {
      assay <- object@assays[[assay_name]]
      
      if (length(assay@layers) == 0) next
      
      valid_layers <- list()
      
      for (lname in names(assay@layers)) {
        layer <- assay@layers[[lname]]
        if (inherits(layer, c("matrix", "dgCMatrix"))) {
          valid_layers[[lname]] <- layer
        }
      }
      
      for (lname in names(assay@layers)) {
        if (!lname %in% names(valid_layers)) {
          warning(sprintf("Removing incompatible layer '%s' from assay '%s'",
                          lname, assay_name), call. = FALSE)
          object@assays[[assay_name]]@layers[[lname]] <- NULL
        }
      }
    }
    
    result <- safe_subset(object)
    
    if (is.null(result)) {
      stop("Subsetting failed even after layer repair")
    }
  }
  
  # Optional preprocessing
  if (!preprocess) {
    return(result)
  }
  
  result <- NormalizeData(result)
  result <- FindVariableFeatures(result)
  
  if (length(unique(Idents(result))) > 1) {
    result <- ScaleData(result)
  } else {
    warning("Skipping ScaleData: only one identity level present", call. = FALSE)
  }
  
  # Dimensional reduction
  if (use.harmony &&
      harmony.group %in% colnames(result@meta.data) &&
      length(unique(result@meta.data[[harmony.group]])) > 1) {
    
    result <- harmony::RunHarmony(result, group.by.vars = harmony.group)
    reduction <- "harmony"
    
  } else {
    result <- RunPCA(result, npcs = ndims)
    reduction <- "pca"
  }
  
  # Graph, clustering
  if ("integrated" %in% names(result@assays)) {
    DefaultAssay(result) <- "integrated"
  } else {
    message("No 'integrated' assay found, using default assay: ", DefaultAssay(result))
  }
  
  result <- FindNeighbors(result, dims = 1:ndims, reduction = reduction)
  result <- FindClusters(result, resolution = resolution)
  
  # Optional UMAP
  if (run.umap) {
    set.seed(123)
    result <- RunUMAP(
      result,
      dims = 1:ndims,
      reduction = reduction,
      min.dist = min.dist,
      spread = spread
    )
  }
  
  DefaultAssay(result) <- "RNA"
  
  return(result)
}


spatial_cellfeat <- function(
    seurat_object,
    features,
    images = c("slice3A"),
    image.alpha = 0,
    pt.size.factor = 50,
    interactive = FALSE,
    alpha = 1,
    ncol = NULL,
    ...
) {
  
  require(ggplot2)  
  require(Seurat)
  
  # Create plot list
  plot_list <- SpatialFeaturePlot(object = seurat_object,
                                  features = features,
                                  interactive = interactive, 
                                  pt.size.factor = pt.size.factor, 
                                  images = images, 
                                  image.alpha = image.alpha, 
                                  alpha = alpha, 
                                  ncol = ncol,
                                  ...) & theme(legend.position = "right",
                                               legend.background = element_blank(),
                                               legend.key.height = unit(0.6, 'cm'),
                                               legend.key.width = unit(0.35, 'cm'),)
  
  return(plot_list)
}




Plot_GO <- function(data,
                    fold_change, 
                    top_n = 10, 
                    specific_go_terms = NULL, 
                    plot_type = "dotplot",
                    theme_type = "minimal",
                    cluster_genes = FALSE,
                    diagonal_order = FALSE,
                    plot_title = "GO Terms Enrichment", 
                    x_title = "Genes", 
                    y_title = "GO Terms",
                    x_angle = 90, 
                    font_size = 10, 
                    wrap_length = 50,
                    alpha_value = 0.8) { 
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(RColorBrewer)
  require(reshape2)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors <- grDevices::colorRampPalette(colors = colors)(100)
  
  data <- data.frame(data)
  data <- na.omit(data)
  data$Description <- stringr::str_to_upper(data$Description)
  
  if (!is.null(specific_go_terms)) {
    data <- data[data$Description %in% specific_go_terms, ]
    if (nrow(data) == 0) stop("No matching GO terms found.")
  } else {
    data <- data %>% arrange(p.adjust) %>% head(top_n)
  }
  
  plot_data <- data %>%
    dplyr::select(Description, geneID, p.adjust) %>%
    dplyr::mutate(Log10_padj = -log10(p.adjust)) %>%
    separate_rows(geneID, sep = "/") %>%
    dplyr::rename(GO_Term = Description, Gene = geneID)
  
  plot_data$GO_Term <- stringr::str_wrap(plot_data$GO_Term, width = wrap_length)
  
  fold_change_df <- data.frame(Gene = names(fold_change), FoldChange = fold_change)
  plot_data <- merge(plot_data, fold_change_df, by = "Gene", all.x = TRUE)
  plot_data <- plot_data %>% filter(!is.na(FoldChange))
  
  if (cluster_genes) {
    gene_order <- plot_data %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarize(order_score = min(as.numeric(factor(GO_Term))), .groups = "drop") %>%
      dplyr::arrange(order_score)
    plot_data$Gene <- factor(plot_data$Gene, levels = gene_order$Gene)
  } else if (diagonal_order) {
    plot_data$GO_Term <- factor(plot_data$GO_Term, levels = rev(unique(plot_data$GO_Term)))
    gene_order_by_go <- plot_data %>% arrange(GO_Term) %>% pull(Gene) %>% unique()
    plot_data$Gene <- factor(plot_data$Gene, levels = gene_order_by_go)
  }
  
  if (plot_type == "dotplot") {
    p <- ggplot(plot_data, aes(x = Gene, y = GO_Term)) +
      geom_point(aes(size = Log10_padj, fill = FoldChange),
                 shape = 21, color = "black", alpha = alpha_value) +
      scale_size(range = c(1, 4), name = "-Log10(padj)") +
      scale_fill_gradientn(
        colors = colors,
        limits = c(min(fold_change, na.rm = TRUE), max(fold_change, na.rm = TRUE)),
        name = "LogFC",
        guide = guide_colorbar(
          barwidth = 0.7,
          barheight = 5,
          frame.colour = "black",
          ticks.colour = "black"
        )
      ) +
      plot_theme(theme_type = theme_type, x_angle = x_angle, font_size = font_size) +
      labs(title = plot_title, x = x_title, y = y_title)
    
  } else if (plot_type == "heatmap") {
    heatmap_matrix <- reshape2::dcast(plot_data, GO_Term ~ Gene, value.var = "FoldChange", fill = NA)
    rownames(heatmap_matrix) <- heatmap_matrix$GO_Term
    heatmap_matrix <- heatmap_matrix[, -1]
    heatmap_long <- reshape2::melt(as.matrix(heatmap_matrix), varnames = c("GO_Term", "Gene"), value.name = "FoldChange")
    
    if (cluster_genes) {
      heatmap_long$Gene <- factor(heatmap_long$Gene, levels = gene_order$Gene)
    } else if (diagonal_order) {
      heatmap_long$GO_Term <- factor(heatmap_long$GO_Term, levels = rev(unique(heatmap_long$GO_Term)))
      gene_order_by_go <- heatmap_long %>% arrange(GO_Term) %>% pull(Gene) %>% unique()
      heatmap_long$Gene <- factor(heatmap_long$Gene, levels = gene_order_by_go)
    }
    
    p <- ggplot(heatmap_long, aes(x = Gene, y = GO_Term, fill = FoldChange)) +
      geom_tile(color = "white", alpha = alpha_value) +
      scale_fill_gradientn(
        colors = colors,
        limits = c(min(fold_change, na.rm = TRUE), max(fold_change, na.rm = TRUE)),
        na.value = "white",
        name = "LogFC",
        guide = guide_colorbar(
          barwidth = 0.7,
          barheight = 5,
          frame.colour = "black",
          ticks.colour = "black"
        )
      ) +
      plot_theme(theme_type = theme_type, x_angle = x_angle, font_size = font_size) +
      labs(title = plot_title, x = x_title, y = y_title)
  } else {
    stop("Invalid plot_type. Choose 'dotplot' or 'heatmap'.")
  }
  
  return(p)
}



#' Summarize Seurat Object Metadata
#'
#' This function summarizes metadata from a Seurat object and generates a bar plot 
#' for visualization. It allows grouping and conditioning by metadata columns, 
#' supports color customization, and provides options for reordering bars.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param group.by A string specifying the metadata column to group by (default: "sample").
#' @param condition.by A string specifying a second metadata column for conditional grouping (default: NULL).
#' @param pal.setup A string specifying the color palette from RColorBrewer (default: "Set1").
#' @param custom.colors A vector of custom colors to use instead of the default palette (default: NULL).
#' @param use.discrete.colors Logical; if TRUE, uses discrete colors from Set1 (default: FALSE).
#' @param angle.x An integer specifying the angle of x-axis labels (default: 90).
#' @param vjust.x Vertical justification of x-axis labels (default: NULL).
#' @param hjust.x Horizontal justification of x-axis labels (default: NULL).
#' @param x.title A string for the x-axis title (default: same as `group.by`).
#' @param remove.axis.x.text Logical; if TRUE, removes x-axis text (default: FALSE).
#' @param reorder.bars Logical; if TRUE, reorders bars based on a specified column (default: FALSE).
#' @param reorder.column A string specifying the column to reorder bars by (default: NULL).
#' @param plot.variable A string specifying the metric to plot: "total_counts", "total_transcripts", 
#'  "mean_transcripts", "median_transcripts", or "proportion" (default: "total_counts").
#' @param title A string specifying the plot title (default: same as `plot.variable`).
#' @param legend Logical; if FALSE, removes the legend (default: TRUE).
#' @param text.size An integer specifying the text size in the plot (default: 10).
#'
#' @return A list containing:
#'   \item{summary_table}{A summary table with calculated statistics.}
#'   \item{plot}{A ggplot2 object representing the bar plot.}
#'
#' @import Seurat dplyr ggplot2 RColorBrewer
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- summarize_seurat(object = object, group.by = "cell_type", plot.variable = "total_counts")
#'   print(result$plot)
#' }
summarize_seurat <- function(object, 
                             group.by = "sample",
                             assay = "RNA",
                             slot = "counts",
                             condition.by = NULL,  
                             pal.setup = "Set1",
                             custom.colors = NULL,
                             use.discrete.colors = FALSE,
                             theme = "classic", 
                             x.angle = 45, 
                             x.title = group.by,
                             remove.axis.x.text = FALSE,
                             reorder.bars = FALSE,
                             reorder.column = NULL,
                             plot.variable = c("total_transcripts", "mean_transcripts", "median_transcripts", 
                                               "num_expressed_genes", "mean_genes_per_cell", "num_cells"),
                             plot.title = plot.variable,
                             legend = TRUE,
                             legend.title = NULL,
                             legend.key.size = 0.3,
                             legend.position = "right",
                             legend.ncol = 1,
                             legend.alpha = 1,
                             text.size = 10,
                             remove.x.labels = FALSE,
                             remove.x.title = FALSE,
                             remove.y.title = FALSE,
                             y.title = NULL) {
  
  # Load required libraries
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(RColorBrewer)
  
  # Validate the plotting method choice
  plot.variable <- match.arg(plot.variable)
  
  # Extract metadata and gene expression matrix
  metadata <- object@meta.data
  gene_counts <- object[[assay]][slot]
  
  # Ensure group.by exists
  if (!group.by %in% colnames(metadata)) {
    stop(paste("The specified group.by column", group.by, "does not exist in the Seurat object's metadata."))
  }
  
  # Ensure condition.by exists if provided
  if (!is.null(condition.by) && !condition.by %in% colnames(metadata)) {
    stop(paste("The specified condition.by column", condition.by, "does not exist in the Seurat object's metadata."))
  }
  
  # Compute expressed genes per cell
  expressed_genes_per_cell <- Matrix::colSums(gene_counts > 0)  # Number of genes expressed per cell
  metadata$num_expressed_genes <- expressed_genes_per_cell
  
  # Group data
  if (!is.null(condition.by)) {
    summary_stats <- metadata %>%
      dplyr::group_by(!!sym(group.by), !!sym(condition.by)) %>%
      dplyr::summarise(
        total_counts = dplyr::n(),  
        total_transcripts = sum(nCount_RNA),  
        mean_transcripts = mean(nCount_RNA),  
        median_transcripts = median(nCount_RNA),  
        num_expressed_genes = sum(num_expressed_genes),  
        mean_genes_per_cell = mean(num_expressed_genes),  
        .groups = 'drop'  
      ) %>%
      dplyr::group_by(!!sym(group.by)) %>%
      dplyr::mutate(proportion = total_counts / sum(total_counts)) %>%  
      dplyr::ungroup()
  } else {
    summary_stats <- metadata %>%
      dplyr::group_by(!!sym(group.by)) %>%
      dplyr::summarise(
        total_counts = dplyr::n(),  
        total_transcripts = sum(nCount_RNA),  
        mean_transcripts = mean(nCount_RNA),  
        median_transcripts = median(nCount_RNA),  
        num_expressed_genes = sum(num_expressed_genes),  
        mean_genes_per_cell = mean(num_expressed_genes),  
        .groups = 'drop'  
      )
  }
  
  # Add the number of cells per group (corrected)
  num_cells_per_group <- metadata %>%
    dplyr::group_by(!!sym(group.by)) %>%
    dplyr::summarise(num_cells = n(), .groups = 'drop')
  summary_stats <- left_join(summary_stats, num_cells_per_group, by = group.by)
  
  # Ensure plot.variable exists
  if (!(plot.variable %in% colnames(summary_stats))) {
    stop(paste("The specified plot.variable", plot.variable, "does not exist in the summary statistics table."))
  }
  
  # Define colors
  unique_groups <- unique(summary_stats[[ifelse(is.null(condition.by), group.by, condition.by)]])
  num_groups <- length(unique_groups)
  
  if (!is.null(custom.colors)) {
    if (length(custom.colors) < num_groups) {
      stop(paste("Insufficient colors in custom.colors. Needed:", num_groups, "Provided:", length(custom.colors)))
    }
    colors_to_use <- custom.colors
  } else if (use.discrete.colors) {
    colors_to_use <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_groups)
  } else {
    colors_to_use <- colorRampPalette(RColorBrewer::brewer.pal(min(9, num_groups), pal.setup))(num_groups)
  }
  
  # Reorder bars globally if required
  if (reorder.bars) {
    if (is.null(reorder.column)) {
      reorder.column <- plot.variable  # Default to reordering by the plot variable
    }
    if (!reorder.column %in% colnames(summary_stats)) {
      stop(paste("The specified reorder.column", reorder.column, "does not exist in the summary statistics table."))
    }
    
    # Ensure reorder.column is evaluated correctly
    summary_stats <- summary_stats %>%
      dplyr::mutate(!!reorder.column := as.numeric(summary_stats[[reorder.column]])) %>%
      dplyr::arrange(dplyr::desc(!!sym(reorder.column))) %>%
      dplyr::mutate(!!sym(group.by) := factor(!!sym(group.by), levels = unique(!!sym(group.by))))
  }
  
  # Create base bar plot
  p <- ggplot(summary_stats, aes(x = !!sym(group.by), y = !!sym(plot.variable))) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +  
    plot_theme(theme = theme, x.angle = x.angle, font.size = text.size) +
    theme(legend.position = legend.position,
          legend.key.size = unit(legend.key.size, "cm")) +
    guides(color = guide_legend(ncol = legend.ncol, alpha = legend.alpha, title = legend.title)) +
    labs(title = plot.title, y = plot.variable, x = x.title, colour = paste(legend.title)) 
  
  # Add fill aesthetic
  if (!is.null(condition.by)) {
    p <- p + aes(fill = !!sym(condition.by)) + scale_fill_manual(values = colors_to_use, name = legend.title)
  } else {
    p <- p + aes(fill = !!sym(group.by)) + scale_fill_manual(values = colors_to_use, name = legend.title)
  }
  
  if (remove.x.title) {
    p <- p + theme(axis.title.x = element_blank())
  }
  
  if (remove.y.title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  
  if (remove.x.labels) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  
  if (!is.null(y.title)) {
    p <- p + ylab(y.title)
  }
  
  # Hide legend if specified
  if (!legend) {
    p <- p & NoLegend()
  }
  
  # Return results
  return(list(summary_table = summary_stats, plot = p))
}


#' Create Sankey Plot from Seurat Object
#'
#' This function generates a Sankey plot using the `networkD3` package from a Seurat object,
#' based on selected variables (metadata columns). It visualizes transitions between different 
#' levels or categories across the selected variables.
#'
#' @param object Seurat object containing the metadata.
#' @param selected_vars A character vector of metadata column names to visualize in the Sankey plot. 
#'                       Must contain at least two variables.
#' @param color_list A named list of colors where names correspond to the identity (e.g., cell type) 
#'                   and values correspond to the colors. If NULL, default colors are used.
#' @param sinksRight A boolean value indicating whether the sinks (final nodes) should be positioned on the right 
#'                   side of the plot. Defaults to FALSE.
#' @param fontSize An integer specifying the font size of node labels. Defaults to 13.
#' @param nodeWidth An integer specifying the width of the nodes. Defaults to 40.
#' @param nodePadding An integer specifying the padding between nodes. Defaults to 20.
#'
#' @return A Sankey plot object created using the `networkD3` package, which can be rendered in R Markdown or Jupyter Notebooks.
#'
#' @examples
#' \dontrun{
#' sankey_plot(seurat_obj, selected_vars = c("ann_level_1", "ann_level_2"))
#' }
#'
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import networkD3
#' @import jsonlite
#' @export
sankey_plot <- function(seurat_obj, 
                        selected_vars = c("ann_level_1", "ann_level_2", "ann_level_3"),
                        plot.title = NULL,
                        custom_colors = list(),  
                        text.size = 14,
                        show_counts = TRUE,
                        show_percentages = TRUE,
                        show_labels = TRUE,
                        label.justify = "left", 
                        label.nudge = 0.1,
                        flow.alpha = 0.5,
                        x.label_position = "top",
                        custom_x.labels = NULL,  
                        show_x_axis = TRUE,
                        prevent_overlap = TRUE) {    
  
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(ggsankey)
  require(viridis)  
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Ensure selected variables exist
  missing_cols <- setdiff(selected_vars, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in metadata:", paste(missing_cols, collapse = ", ")))
  }
  
  # Convert data to long format
  df <- metadata %>%
    dplyr::select(all_of(selected_vars)) %>%
    make_long(!!!syms(selected_vars))  
  
  # Calculate total per selected variable level
  dagg <- df %>%
    group_by(x, node) %>%
    tally() %>%
    dplyr::mutate(pct = n / sum(n))  
  
  # Merge back into the main dataframe
  df <- left_join(df, dagg, by = c("x", "node"))
  
  # Unique ID for each node
  df$full_node_name <- paste(df$x, df$node, sep = "_")  
  
  # Build color map
  color_map <- c()
  for (var in selected_vars) {
    if (!is.null(custom_colors[[var]])) {
      color_map <- c(color_map, setNames(custom_colors[[var]], paste(var, names(custom_colors[[var]]), sep = "_")))
    }
  }
  
  # # Default to viridis if no custom colors
  # if (length(color_map) == 0) {
  #   viridis_colors <- viridis::viridis(length(unique(df$full_node_name)))
  #   color_map <- setNames(viridis_colors, unique(df$full_node_name))
  # }
  
  # Default to custom_palette if no custom colors
  if (length(color_map) == 0) {
    n_nodes <- length(unique(df$full_node_name))
    palette_colors <- custom_palette(n_nodes, preset = "base")  # choose any preset you like
    color_map <- setNames(palette_colors, unique(df$full_node_name))
  }
  
  # Check missing nodes
  unique_nodes <- unique(df$full_node_name)
  missing_nodes <- setdiff(unique_nodes, names(color_map))
  if (length(missing_nodes) > 0) {
    stop(paste("Error: Missing colors for nodes:", paste(missing_nodes, collapse = ", ")))
  }
  
  # Assign colors
  color_scale <- scale_fill_manual(values = color_map)
  
  #- ADAPTIVE LABEL COLORS-
  get_text_color <- function(hex_color) {
    rgb_vals <- col2rgb(hex_color) / 255
    brightness <- sum(rgb_vals * c(0.299, 0.587, 0.114))  
    if (brightness > 0.5) "black" else "white"
  }
  df$label_color <- sapply(df$full_node_name, function(name) get_text_color(color_map[name]))
  
  #- Labels-
  df$label <- df$node  
  if (show_labels) {
    if (show_counts & show_percentages) {
      df$label <- paste0(df$node, " n=", df$n, " (", round(df$pct * 100, 2), "%)")
    } else if (show_counts) {
      df$label <- paste0(df$node, " n=", df$n)
    } else if (show_percentages) {
      df$label <- paste0(df$node, " (", round(df$pct * 100, 2), "%)")
    }
  } else {
    df$label <- ""  
  }
  
  # Adjust justification
  hjust_value <- ifelse(label.justify == "left", 1, ifelse(label.justify == "center", 0.5, 0))
  nudge_x_value <- ifelse(label.justify == "left", label.nudge, ifelse(label.justify == "right", -label.nudge, 1))
  
  # Prevent overlap by ordering
  if (prevent_overlap) {
    df <- df %>%
      group_by(x) %>%
      arrange(desc(n)) %>%
      mutate(node = factor(node, levels = unique(node)))
  }
  
  #- Plot-
  pl <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
                       fill = full_node_name, label = label)) +
    geom_sankey(flow.alpha = flow.alpha, node.color = "black", show.legend = TRUE) +
    geom_sankey_label(aes(color = I(label_color)), size = 3, 
                      hjust = hjust_value, nudge_x = nudge_x_value) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(1, 1, 1, 1, "cm"),
          axis.text.x = element_text(color = "black", size = text.size, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_x_discrete(position = x.label_position) +
    color_scale +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.6))) +
    labs(title = plot.title, fill = "Nodes")
  
  # Custom x-axis labels
  if (!is.null(custom_x.labels) & length(custom_x.labels) == length(selected_vars)) {
    pl <- pl + scale_x_discrete(labels = custom_x.labels)
  }
  
  # Hide x-axis if needed
  if (!show_x_axis) {
    pl <- pl + theme(axis.text.x = element_blank())
  }
  
  return(pl)
}



#' Cell map visualization for Seurat objects
#'
#' @param object Seurat object
#' @param group.by Column(s) in metadata to group cells
#' @param reduction Dimensionality reduction to use (default "umap")
#' @param dims Dimensions to plot (1,2) or 3 for 3D
#' @param shuffle Logical, shuffle cells
#' @param raster Logical, rasterize plot
#' @param alpha Point transparency
#' @param repel Logical, repel labels
#' @param n.cells Logical, show number of cells in labels
#' @param label Logical, show cluster labels
#' @param label.size Cluster label size
#' @param label.face Cluster label font face
#' @param colors Named vector of colors
#' @param figplot Logical, minimal figure plot for figure panels
#' @param no.axes Logical, hide axes
#' @param plot.ttl Plot title
#' @param legend Logical, show legend
#' @param leg.ttl Legend title
#' @param item.size Legend item size
#' @param leg.pos Legend position
#' @param leg.just Legend justification
#' @param leg.dir Legend direction
#' @param leg.size Legend font size
#' @param leg.ncol Legend number of columns
#' @param item.border Logical, border around legend items
#' @param font.size Base font size
#' @param pt.size Point size
#' @param dark Logical, dark theme
#' @param total.cells Logical, include total cells in title
#' @param threeD Logical, 3D plot
#' @param theme Theme name
#' @param facet.bg Logical, add facet background
#' @param ... Additional arguments passed to DimPlot
#' @return ggplot or plotly object
#' @export
#' 
cellmap <- function(
    object, 
    group.by = NULL, 
    reduction = "umap", 
    dims = c(1,2), 
    shuffle = FALSE,
    raster = NULL,
    stroke.size = NULL,
    raster.dpi = c(512, 512), 
    alpha = 1, 
    repel = FALSE, 
    n.cells = TRUE,
    label = FALSE, 
    label.size = 3.5, 
    label.face = "plain", 
    cols = NULL, 
    figplot = FALSE, 
    no.axes = FALSE, 
    plot.ttl = NULL, 
    legend = TRUE,
    leg.ttl = NULL, 
    leg.ttl.size = font.size,
    item.size = 3.5, 
    leg.pos = "right", 
    leg.just = "center",
    leg.dir = "vertical", 
    leg.size = 10, 
    leg.ncol = NULL, 
    item.border = TRUE,
    font.size = 10, 
    pt.size = NULL, 
    dark = FALSE, 
    total.cells = FALSE,
    threeD = FALSE, 
    theme.style = "classic", 
    facet.bg = FALSE,
    ...
) {
  
  if (!is.null(list(...)$theme)) {
    theme.style <- list(...)$theme
  }
  
  # Helper: prepare object and colors
  .prepare_object <- function(obj, group, cols=NULL){
    stopifnot(inherits(obj,"Seurat"))
    if(is.null(group)) group <- "ident"
    
    if(group=="ident"){
      obj@meta.data$ident <- Idents(obj)
    } else if(!group %in% colnames(obj@meta.data)){
      stop(paste("Grouping column", group, "not found."))
    }
    
    values <- as.character(obj@meta.data[[group]])
    values[is.na(values)] <- "Unknown"
    levels_group <- if(is.factor(obj@meta.data[[group]])) levels(obj@meta.data[[group]]) else unique(values)
    if("Unknown" %in% values & !("Unknown" %in% levels_group)) levels_group <- c(levels_group, "Unknown")
    obj@meta.data[[group]] <- factor(values, levels=levels_group)
    Idents(obj) <- obj@meta.data[[group]]
    
    # # Drop unused levels
    # obj@meta.data[[group]] <- droplevels(factor(obj@meta.data[[group]]))
    # values <- as.character(obj@meta.data[[group]])
    # values[is.na(values)] <- "Unknown"
    # obj@meta.data[[group]] <- factor(values)
    # Idents(obj) <- obj@meta.data[[group]]
    # levels_group <- levels(obj@meta.data[[group]])
    
    if(is.null(cols)){
      cols <- custom_palette(length(levels_group))
      names(cols) <- levels_group
      if("Unknown" %in% levels_group) cols["Unknown"] <- "gray70"
    } else {
      if(is.null(names(cols))) cols <- setNames(cols[seq_along(levels_group)], levels_group)
      missing <- setdiff(levels_group, names(cols))
      if(length(missing)) cols[missing] <- "gray70"
      cols <- cols[levels_group]
    }
    
    list(obj=obj, cols=cols, levels_group=levels_group)
  }
  
  # 3D plotting helper
  .plot_3D <- function(obj, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells){
    emb <- obj@reductions[[reduction]]@cell.embeddings
    df <- data.frame(x=emb[, dims[1]], y=emb[, dims[2]], z=emb[, dims[3]], cluster=Idents(obj))
    hover_labels <- if(n.cells){
      tbl <- table(df$cluster)
      paste0(df$cluster, " (", tbl[as.character(df$cluster)], ")")
    } else as.character(df$cluster)
    
    p3d <- plotly::plot_ly(df, x=~x, y=~y, z=~z, color=~cluster, colors=cols,
                           type="scatter3d", mode="markers",
                           marker=list(size=pt.size, opacity=alpha, line=list(width=0)),
                           text=hover_labels, hoverinfo="text")
    if(label){
      centers <- df %>% dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), z=median(z))
      p3d <- p3d %>% plotly::add_text(data=centers, x=~x, y=~y, z=~z, text=~cluster, textposition="top center")
    }
    p3d
  }
  
  if(length(group.by) > 1){
    plots <- lapply(group.by, function(g){
      cellmap(object = object, group.by = g, shuffle = shuffle,raster = raster,stroke.size = stroke.size,alpha = alpha,
              repel = repel,reduction = reduction,dims = dims,n.cells = n.cells,label = label,label.size = label.size,
              label.face = label.face,cols = cols,figplot = figplot,plot.ttl = g,legend = legend,leg.ttl = g,item.size = item.size,
              leg.pos = leg.pos,leg.just = leg.just,leg.dir = leg.dir,leg.ncol = leg.ncol,font.size = font.size,item.border = item.border,
              pt.size = pt.size,dark = dark,total.cells = total.cells,threeD = threeD,theme.style = theme.style,facet.bg = facet.bg,
              ...
      )
    })
    return(patchwork::wrap_plots(plots))
  }
  # Prepare object & colors
  prep <- .prepare_object(object, group.by, cols)
  object <- prep$obj
  cols <- prep$cols
  levels_group <- prep$levels_group
  if(is.null(leg.ncol)) leg.ncol <- if(length(levels_group) > 18) 2 else 1
  
  # Validate reduction
  if(!(reduction %in% names(object@reductions))){
    stop(paste0("Reduction '", reduction, "' not found. Available: ", paste(names(object@reductions), collapse=", ")))
  }
  emb <- object@reductions[[reduction]]@cell.embeddings
  if(max(dims) > ncol(emb)) stop("Selected dims exceed available dimensions in reduction.")
  
  # 3D plotting
  if(threeD || length(dims) == 3) return(.plot_3D(object, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells))
  
  # 2D plotting
  plt <- Seurat::DimPlot(object, group.by=group.by, shuffle=shuffle, raster=raster, pt.size=pt.size,
                         repel=repel, alpha=alpha, reduction=reduction, dims=dims, raster.dpi=raster.dpi,...)
  
  # Auto raster for large datasets
  if(is.null(raster)) raster <- ncol(object) > 50000
  
  # Adjust point size for large datasets
  if(is.null(pt.size)) pt.size <- if(ncol(object) > 50000) 0.3 else 0.5
  
  # reset Seurat's forced theme_classic
  plt <- plt + ggplot2::theme_void()
  
  present_levels <- levels(droplevels(object@active.ident))
  
  if (n.cells) {
    cell.nb <- table(object@active.ident)[present_levels]
    clust.lab <- paste0(present_levels, " (", cell.nb, ")")
  } else {
    clust.lab <- present_levels
  }
  
  cols_use <- cols[present_levels]
  
  leg.ttl <- if(is.null(leg.ttl)) group.by else leg.ttl
  
  plt <- plt + ggplot2::scale_color_manual(
    breaks = present_levels,
    labels = clust.lab,
    values = cols_use
  )
  
  
  if (legend) {
    plt <- plt & ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = if(item.border)
          list(size = item.size, shape = 21, color = "black",
               stroke = 0.2, fill = unname(cols))
        else
          list(size = item.size),
        ncol = leg.ncol,
        title = leg.ttl,
        keyheight = grid::unit(0.25,"cm"),
        keywidth  = grid::unit(0.25,"cm")
      )
    )
  } else {
    plt <- plt & ggplot2::guides(color = "none")
  }
  
  # Plot title with total cells
  if(total.cells){
    plot.ttl <- paste0(plot.ttl, " (n=", format(ncol(object), big.mark=","), ")")
  }
  #if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl)
  if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl) else plt <- plt + labs(title = NULL)
  
  # Set default legend title size
  if (is.null(leg.ttl.size)) leg.ttl.size <- font.size
  
  # Apply plot_theme using do.call
  theme_args <- list(
    theme.style = theme.style,
    font.size = font.size,
    leg.size  = leg.size,
    leg.pos   = leg.pos,
    leg.dir   = leg.dir,
    leg.ttl   = leg.ttl,
    leg.ttl.size = leg.ttl.size,
    facet.bg  = facet.bg,
    mode      = if(dark) "dark" else "light"
  ) 
  
  if(figplot){
    # Warn if the user specified a non-classic theme
    if(!missing(theme.style) && theme.style != "classic"){
      warning(sprintf(
        "figplot = TRUE ignores custom themes (theme.style = '%s'). Use figplot = FALSE for full theming.",
        theme.style
      ))
    }
    
    # Apply plot_theme for figplot figure
    plt <- plt &
      do.call(plot_theme, c(
        theme_args,
        list(
          x.ttl = FALSE,
          ticks = FALSE,
          line = FALSE,
          border = FALSE,
          grid.major = FALSE,
          grid.minor = FALSE,
          panel.fill = "white"
        ),
        list(...)
      ))    
    text_col <- "black"
  } else {
    # Apply plot_theme normally
    plt <- plt & do.call(plot_theme, c(theme_args, list(...)))
  }
  
  
  # Add cluster labels if requested
  if(label){
    umap_data <- dplyr::tibble(x=emb[, dims[1]], y=emb[, dims[2]], cluster=as.character(object@active.ident)) %>%
      dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups="drop")
    plt <- plt + ggrepel::geom_text_repel(
      data=umap_data, aes(x, y, label=cluster),
      color = if(dark) "white" else "black",
      fontface = label.face,
      bg.color = if(dark) "#3A3A3A" else "grey95",
      bg.r = 0.1, size = label.size, seed = 42
    )
  }
  
  # figplot arrow axes (minimal figure)
  if(figplot){
    x.lab.reduc <- plt$labels$x %||% paste0(toupper(reduction), dims[1])
    y.lab.reduc <- plt$labels$y %||% paste0(toupper(reduction), dims[2])
    plt <- plt & Seurat::NoAxes()
    L <- 0.12
    axis.df <- data.frame(x0=c(0,0), y0=c(0,0), x1=c(L,0), y1=c(0,L))
    axis.plot <- ggplot2::ggplot(axis.df) +
      ggplot2::geom_segment(ggplot2::aes(x=x0, y=y0, xend=x1, yend=y1), linewidth=0.4, lineend="round") +
      ggplot2::xlab(x.lab.reduc) + ggplot2::ylab(y.lab.reduc) +
      ggplot2::coord_fixed() + ggplot2::theme_classic(base_size=font.size) +
      ggplot2::theme(plot.background=ggplot2::element_rect(fill="transparent", colour=NA),
                     panel.background=ggplot2::element_rect(fill="transparent", colour=NA),
                     axis.text=ggplot2::element_blank(),
                     axis.ticks=ggplot2::element_blank(),
                     axis.line=ggplot2::element_blank(),
                     panel.border=ggplot2::element_blank(),
                     axis.title=ggplot2::element_text(size=font.size, face="plain"),
                     plot.margin=ggplot2::margin(0,0,0,0))
    figure.layout <- c(patchwork::area(t=1,l=1,b=11,r=11), patchwork::area(t=10,l=1,b=11,r=2))
    return(plt + axis.plot + patchwork::plot_layout(design=figure.layout))
  }
  
  
  if(!legend) plt <- plt & Seurat::NoLegend()
  if(no.axes) plt <- plt & Seurat::NoAxes()
  
  plt
}

#' Cell Dot Plot (Enhanced Seurat DotPlot)
#'
#' A cleaner, more customizable wrapper around **Seurat::DotPlot**, providing
#' improved color handling, optional dot outlines, flexible axis formatting,
#' legend placement, and theme control. Useful for visualizing gene expression
#' patterns across clusters or metadata-defined groups.
#'
#' @param object A Seurat object.
#' @param features Character vector of features (genes or metadata fields) to plot.
#' @param group.by Column in `object@meta.data` used to group cells.
#'   Default: `"seurat_clusters"`.
#'
#' @param th.cols Color palette name from **RColorBrewer** used for the
#'   expression gradient. Default: `"Reds"`.
#' @param rev.th.cols Logical; reverse the gradient palette. Default: FALSE.
#'
#' @param dot.scale Numeric scale factor controlling the dot size range.
#'   Passed to `Seurat::DotPlot`. Default: 4.5.
#' @param dot.outline Logical; draw outlines around dots. Default: FALSE.
#'
#' @param x.angle Angle for x-axis labels (degrees). Default: 90.
#' @param vjust.x,hjust.x Vertical and horizontal justification for x labels.
#'
#' @param flip Logical; swap x and y axes using `coord_flip()`. Default: FALSE.
#'
#' @param font.size Base font size passed to internal theme helper. Default: 8.
#'
#' @param plot.title Optional plot title.
#'
#' @param leg.size Legend text size. Default: 8.
#' @param leg.pos Position of the legend (e.g., `"right"`, `"bottom"`).
#' @param leg.just Legend justification.
#' @param leg.hjust Logical; if TRUE, use a horizontal legend layout when possible.
#'
#' @param x.axis.pos Position of the x-axis (`"top"` or `"bottom"`).
#' @param theme ggplot2 theme name used by the internal theme helper.
#'
#' @param x.face,y.face Logical; italic styling for x and/or y-axis labels.
#' @param x.ttl,y.ttl Logical; italic styling for x and/or y-axis titles.
#'
#' @param ... Additional parameters passed to `Seurat::DotPlot()`.
#'
#' @details
#' This function enhances the standard Seurat dot plot by providing:
#' * Customizable Brewer color gradients  
#' * Optional dot outlines  
#' * Flexible axis label styling  
#' * Improved legend customization and ordering  
#' * Optional axis flipping  
#'
#' It retains all functionality of `Seurat::DotPlot` while adding cleaner,
#' publication-ready defaults.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' celldot(pbmc, features = c("MS4A1","CD3D"))
#'
#' celldot(pbmc, features = c("MS4A1","CD14"), th.cols = "Blues",
#'         dot.outline = TRUE, flip = TRUE)
#' }
#'
#' @export
#' 
celldot <- function(object, features, group.by="seurat_clusters", th.cols="Reds",
                    rev.th.cols=FALSE, dot.scale=4.5, x.angle=90, vjust.x=NULL,
                    hjust.x=NULL, flip=FALSE, font.size=8, plot.title=NULL,
                    leg.size=10, leg.pos="right", leg.just="bottom", leg.hjust=FALSE,
                    x.axis.pos="bottom", theme="classic", x.face=FALSE, y.face=FALSE,
                    x.ttl=FALSE, y.ttl=FALSE, dot.outline=FALSE, ...) {
  
  stopifnot(inherits(object,"Seurat"))
  object <- Seurat::SetIdent(object, value=group.by)
  features <- unique(features)
  
  pal <- RColorBrewer::brewer.pal(9, th.cols)
  if (rev.th.cols) pal <- rev(pal)
  
  outline_col <- if (dot.outline) "gray60" else NA
  outline_stroke <- if (dot.outline) 0.5 else 0
  
  plt <- suppressWarnings({
    suppressMessages({
      Seurat::DotPlot(object, features=features, dot.scale=dot.scale, ...) 
    })
  }) +
    scale_color_gradientn(colors=pal, oob=scales::squish) +
    geom_point(aes(size=pct.exp), shape=21, colour=outline_col, stroke=outline_stroke) +
    labs(title=plot.title, color="Average\nExpression", size="Percent\nExpressed") +
    plot_theme(theme=theme, font.size=font.size, x.angle=x.angle,
               x.hjust=hjust.x, x.vjust=vjust.x, xy.val=TRUE, x.lab=TRUE, y.lab=TRUE, ...) +
    theme(
      axis.text.x = if (x.face || (flip && y.face)) element_text(face="italic") else element_text(),
      axis.text.y = if (y.face || (flip && x.face)) element_text(face="italic") else element_text(),
      axis.title = element_blank(), 
      legend.spacing.y = unit(0.05, "cm"),
      legend.spacing.x = unit(0.05, "cm"),
      legend.box.spacing = unit(0.05, "cm"),
      legend.margin = margin(2,2,2,2)
    )
  
  if (flip) plt <- plt + coord_flip()
  if (is.list(features)) plt <- plt + theme(strip.text.x=element_text(angle=45))
  
  # Legend positioning
  # Legend positioning outside plot, bottom-right
  if (!is.null(leg.pos)) {
    if (leg.pos == "right") {
      plt <- plt + theme(
        legend.position = "right",
        legend.justification = c("right","bottom"),
        legend.box.just = "right",
        legend.box.margin = margin(0,0,0,0)
      )
    } else if (leg.pos == "left") {
      plt <- plt + theme(
        legend.position = "left",
        legend.justification = c("left","bottom"),
        legend.box.just = "left",
        legend.box.margin = margin(0,0,0,0)
      )
    } else if (leg.pos == "top") {
      plt <- plt + theme(
        legend.position = "top",
        legend.justification = c("right","top"),
        legend.box.just = "right"
      )
    } else if (leg.pos == "bottom") {
      plt <- plt + theme(
        legend.position = "bottom",
        legend.justification = c("right","bottom"),
        legend.box.just = "right"
      )
    }
  }
  
  # Keep original guides for color and size
  guide_color <- guide_colorbar(frame.colour="black", ticks.colour="black")
  guide_size <- guide_legend(override.aes=list(shape=21, colour=outline_col, fill="black"))
  guide_color$order <- 1
  guide_size$order  <- 2
  plt <- plt + guides(color=guide_color, size=guide_size)
  
  plt
}



#' Cell Feature Violin Plot with Statistics
#'
#' Plots expression or metadata features as violin plots for a Seurat object,
#' optionally adding median points, shared y-axis scaling, flipped axes, and
#' statistical comparisons (Wilcoxon for 2 groups, Kruskal-Wallis for >2 groups).
#'
#' @param obj A Seurat object.
#' @param features Character vector of feature names (genes or metadata columns) to plot.
#' @param ncol Number of columns in the output patchwork plot. Defaults to sqrt(#features / 1.5).
#' @param stack Logical; if TRUE, plots are stacked in a single column.
#' @param shared.y Logical; if TRUE, all violins share the same y-axis.
#' @param ttl.pos Position of subplot titles: "center", "left", or "right".
#' @param group.by Metadata column to group by. Defaults to "seurat_clusters".
#' @param split.by Optional metadata column to split violins by.
#' @param assay Assay to pull data from. Default is "RNA".
#' @param slot Slot to use for expression values. One of "data", "counts", or "scale.data".
#' @param log Logical; if TRUE, log-transform the expression values.
#' @param cols Optional named vector of colors for each group. If NULL, defaults are used.
#' @param med Logical; if TRUE, overlay median points on each violin.
#' @param med.size Size of median points if med = TRUE.
#' @param pt.size Size of jittered points. Set to 0 to hide points.
#' @param border.size Size of the violin border lines.
#' @param font.size Base font size for titles and labels.
#' @param theme ggplot2 theme to use: "classic", "minimal", etc.
#' @param x.angle Rotation angle of x-axis labels.
#' @param leg.pos Position of legend: "none", "right", "left", etc.
#' @param title Optional overall title for the patchwork plot.
#' @param rm.subtitles Logical; if TRUE, removes individual subplot titles.
#' @param flip Logical; if TRUE, flips x and y axes.
#' @param auto.resize Logical; if TRUE, sets dynamic width/height attributes.
#' @param ylab.global Global y-axis label. Defaults to expression level.
#' @param xlab.global Global x-axis label. Defaults to blank.
#' @param pairwise Logical; if TRUE, perform pairwise comparisons between groups.
#' @param add.stats Logical; if TRUE, add p-values to plots.
#' @param show.pval Logical; if TRUE, show p-values above violins.
#' @param pval.label Character; label type for p-values, e.g., "p.signif" or "p.format".
#' @param ... Additional arguments passed to ggplot2 layers.
#'
#' @return A patchwork object containing the violin plots.
#' @examples
#' \dontrun{
#' cellvio(sub, features = c("MEG3","TP63","HES6"),
#'         group.by = "ann_level_2",
#'         pt.size = 0.1, ncol = 3, pairwise = TRUE,
#'         font.size = 10, show.pval = TRUE)
#' }
#' @export
#' 
cellvio <- function(
    obj, features,
    ncol = NULL,
    shared.y = FALSE,
    ttl.pos = c("center", "left", "right"),
    group.by = "seurat_clusters",
    split.by = NULL,
    stack = FALSE,
    assay = "RNA",
    slot = "data",
    log = FALSE,
    cols = NULL,
    med = FALSE,
    med.size = 1,
    pt.size = 0,
    border.size = 0.1,
    theme = "classic",
    leg.pos = "none",
    x.angle = 45,
    title = NULL,
    rm.subttl = FALSE,
    flip = FALSE,
    auto.resize = TRUE,
    ylab.global = NULL,
    xlab.global = NULL,
    add.stats = FALSE,
    show.pval = FALSE,
    pairwise = FALSE,
    pval.label = "p.signif",
    font.size = 10,
    ...
) {
  
  stopifnot(inherits(obj, "Seurat"))
  if (length(features) == 0) stop("features must be provided.")
  ttl.pos <- match.arg(ttl.pos)
  
  
  # determine ncol
  if (is.null(ncol)) {
    ncol <- if (stack) 1 else max(1, ceiling(sqrt(length(features) / 1.5)))
  }
  
  # Set identities if grouping
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(obj@meta.data))
      stop(paste(group.by, "not found in metadata."))
    Idents(obj) <- group.by
  }
  
  # Determine which features exist
  f.expr <- intersect(features, rownames(obj[[assay]]))
  f.meta <- intersect(features, colnames(obj@meta.data))
  features <- unique(c(f.expr, f.meta))
  if (!length(features)) stop("No features found in assay or metadata.")
  
  # Shared y-scale
  ymax <- NULL
  if (shared.y) {
    vals <- c(
      if (length(f.expr)) as.numeric(Seurat::GetAssayData(obj, assay, leyer)[f.expr, ]),
      if (length(f.meta)) as.numeric(as.matrix(obj@meta.data[, f.meta, drop = FALSE]))
    )
    vals <- vals[is.finite(vals)]
    if (length(vals)) ymax <- max(vals)
  }
  
  # Colors
  if (is.null(cols)) {
    g <- tryCatch(unique(obj[[group.by]][, 1]), error = \(e) NULL)
    cols <- custom_palette(length(g))
    #cols <- scales::hue_pal()(if (is.null(g)) 8 else length(g))
    #cols <- ggpubr::get_palette("npg", length(g))
    names(cols) <- g
  }
  
  # Helper to build a single violin with stats
  vln <- function(f) {
    
    if (log && slot == "counts") {
      obj[[assay]]@data[f, ] <- log1p(obj[[assay]]@counts[f, ])
      slot_use <- "data"
    } else {
      slot_use <- slot
    }
    
    df <- Seurat::FetchData(obj, vars = c(group.by, f))
    names(df)[2] <- "value"
    df[[group.by]] <- factor(df[[group.by]]) # ensure factor
    
    # Base plot
    # p <- ggplot(df, aes_string(group.by, "value", fill = group.by)) +
    #   geom_violin(scale = "width", color = "black", size = border.size) +
    #   scale_fill_manual(values = cols)
    
    p <- suppressWarnings({
      suppressMessages({Seurat::VlnPlot(
        obj,
        features = f,
        group.by = group.by,
        split.by = split.by,
        assay = assay,
        slot = slot,
        pt.size = pt.size,
        cols = cols,
        ...
      ) + scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.25))
      )
      })
    })
    
    # Points
    #if (pt.size > 0) p <- p + geom_jitter(width = 0.1, size = pt.size, alpha = 0.6)
    
    # Add statistics
    if (add.stats && show.pval && length(levels(df[[group.by]])) > 1) {
      cmp <- if (pairwise) utils::combn(levels(df[[group.by]]), 2, simplify = FALSE) else NULL
      stat_fun <- if (length(levels(df[[group.by]])) == 2) "wilcox.test" else "kruskal.test"
      y_max <- max(df$value, na.rm = TRUE)
      y_step <- (ymax %||% y_max) * 0.05
      p <- p + ggpubr::stat_compare_means(
        method = stat_fun,
        comparisons = cmp,
        label = pval.label,
        hide.ns = FALSE,
        label.y = y_max + seq(0, by = y_step, length.out = length(cmp))
      )
    }
    
    
    # Titles
    p <- if (!rm.subttl) p + labs(title = f) else p + labs(title = NULL)
    
    # Make feature titles bold + italic
    p <- p + theme(
      plot.title = element_text(face = "bold.italic")
    )
    
    .style_layers <- function(
    p,
    violin_lw = 0.15,
    point_size = NULL,
    jitter_width = NULL
    ) {
      for (i in seq_along(p$layers)) {
        
        layer <- p$layers[[i]]
        
        # Violin outline 
        if (inherits(layer$geom, "GeomViolin")) {
          layer$aes_params$linewidth <- violin_lw
        }
        
        # Points (Seurat uses GeomPoint + position_jitterdodge) 
        if (inherits(layer$geom, "GeomPoint")) {
          
          if (!is.null(point_size)) {
            layer$aes_params$size  <- point_size
            layer$aes_params$alpha <- 0.6
          }
          
          if (!is.null(jitter_width) &&
              inherits(layer$position, "PositionJitterdodge")) {
            layer$position$width <- jitter_width
          }
        }
        
        p$layers[[i]] <- layer
      }
      p
    }
    
    p <- .style_layers(
      p,
      violin_lw  = border.size,
      point_size = if (pt.size > 0) pt.size else NULL,
      jitter_width = 0.08
    )
    
    # Theme & formatting
    p <- p + plot_theme(theme = theme, font.size = font.size, x.angle = x.angle,
                        leg.pos = leg.pos, x.ttl = FALSE, ttl.pos = ttl.pos,...) +
      theme(
        plot.title = element_text(face = "bold.italic")
      )
    
    # Force final legend position
    if (!is.null(leg.pos)) {
      p <- p + theme(legend.position = leg.pos)
    }
    
    if (med) p <- p + stat_summary(fun = median, geom = "point", shape = 3, size = med.size)
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    if (flip) p <- p + coord_flip()
    p + ylab(NULL)
  }
  
  # number of cols
  if (is.null(ncol))
    ncol <- max(1, ceiling(sqrt(length(features) / 1.5)))
  
  plist <- lapply(features, vln)
  total <- length(plist)
  
  # Only show x-axis on bottom plots
  bottom <- sapply(1:ncol, \(i) max(seq(i, total, by = ncol)))
  for (i in seq_along(plist)) {
    if (!(i %in% bottom)) {
      plist[[i]] <- plist[[i]] +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }
  }
  
  # Layout
  combo <- patchwork::wrap_plots(plist, ncol = ncol) +
    patchwork::plot_layout(guides = "collect")
  
  if (!is.null(title)) {
    combo <- combo +
      patchwork::plot_annotation(
        title = title,
        theme = theme(title = element_text(face = "bold"))
      )
  }
  
  # Global labels
  auto_y <- switch(slot,
                   data = "Expression level",
                   counts = "Raw counts",
                   scale.data = "Scaled expression",
                   "Expression level")
  
  ylab <- if (is.null(ylab.global)) auto_y else ylab.global
  xlab <- if (is.null(xlab.global)) "" else xlab.global
  
  plt <- cowplot::ggdraw(combo) +
    cowplot::draw_label(ylab, x = -0.01, y = 0.55, angle = 90, size = font.size) +
    cowplot::draw_label(xlab, x = 0.5, y = 0.02, size = font.size) +
    theme(plot.margin = margin(15, 15, 15, 15))
  
  # Auto resize attributes
  if (auto.resize) {
    ng <- length(unique(obj[[group.by]][, 1]))
    attr(plt, "dynamic_width") <- 6 + ng * 0.3
    attr(plt, "dynamic_height") <- 4 + length(features) * 0.25
  }
  
  plt
}

#' Feature Plot Wrapper for Seurat Objects
#'
#' This function creates feature plots for one or multiple features in a Seurat object.
#' It preserves the default Seurat behavior while providing options for custom color palettes,
#' viridis palettes, NA color, legend merging, orientation, and size.
#'
#' @param seurat_object A Seurat object.
#' @param features A character vector of feature names to plot (genes or metadata columns).
#' @param colors_use Optional vector of colors to use for the gradient. Must have at least 2 colors.
#' @param theme_color Color theme to use if `colors_use` is NULL. Default is "Reds".
#' @param use_viridis Logical, whether to use a viridis palette. Default is FALSE.
#' @param viridis_option Character, the viridis palette option to use ("viridis", "magma", "plasma", "inferno", "cividis", etc.). Default is "viridis".
#' @param reverse_colors Logical, whether to reverse the color scale. Default is FALSE.
#' @param na_color Color to use for NA values. Default is "lightgray".
#' @param na_cutoff Numeric cutoff for minimal expression to show. Default is 1e-9.
#' @param order Logical, whether to order points by expression. Default is FALSE.
#' @param pt.size Numeric, size of points. Default automatically adjusts to number of cells.
#' @param base_size Numeric, base font size for plot title and axes. Default is 14.
#' @param legend.text.size Numeric, font size for legend text. Default is 12.
#' @param reduction Character, which dimensional reduction to use. Default uses Seurat's `DefaultDimReduc()`.
#' @param raster Logical, whether to rasterize the plot for speed. Default TRUE for >2e5 cells.
#' @param raster.dpi Numeric vector of length 2, resolution for rasterized plots. Default c(512,512).
#' @param split.by Character, metadata column to split plots by.
#' @param ncol Number of columns when combining plots. Default NULL.
#' @param layer Character, which assay slot to use. Default is "data".
#' @param label Logical, whether to add cluster labels. Default FALSE.
#' @param no_axes Logical, whether to remove axes. Default FALSE.
#' @param blend Logical, whether to create a blend plot for exactly 2 features. Default FALSE.
#' @param merge_legend Logical, whether to merge legends for multiple plots. Default FALSE.
#' @param merge.leg.pos Position of merged legend ("right", "bottom", etc.). Default "right".
#' @param legend.orientation Legend orientation: "vertical" or "horizontal". Default "vertical".
#' @param legend.width Unit object specifying width of legend bar. Default `unit(1,"mm")`.
#' @param legend.height Unit object specifying height of legend bar. Default `unit(8,"mm")`.
#' @param ... Additional arguments passed to Seurat's `FeaturePlot()`.
#'
#' @return A ggplot object (or patchwork) containing the feature plots.
#' @export
#'
#' @examples
#' \dontrun{
#' cellfeat(
#'   seurat_object = obj,
#'   features = c("TEX15","TP63"),
#'   merge_legend = TRUE,
#'   merge.leg.pos = "bottom",
#'   legend.orientation = "horizontal",
#'   legend.width = unit(10,"mm"),
#'   legend.height = unit(1,"mm"),
#'   theme_color = "Reds",
#'   no_axes = TRUE,
#'   ncol = 1,
#'   na_cutoff = 0.5
#' )
#' }
cellfeat <- function(
    seurat_object,
    features,
    colors_use = NULL,
    theme_color = "Reds",
    use_viridis = FALSE,
    viridis_option = "viridis",
    reverse_colors = FALSE,
    na_color = "lightgray",
    na_cutoff = 1e-9,
    order = FALSE,
    pt.size = NULL,
    base_size = 10,
    legend.text.size = 10,
    reduction = NULL,
    raster = NULL,
    raster.dpi = c(512,512),
    split.by = NULL,
    ncol = NULL,
    layer = "data",
    label = FALSE,
    no_axes = FALSE,
    blend = FALSE,
    merge_legend = FALSE,
    merge.leg.pos = "right",
    legend.orientation = c("vertical","horizontal"),
    legend.title.pos = "left",
    legend.title.angle = 90,
    legend.width = grid::unit(0.8,"mm"),
    legend.height = grid::unit(6,"mm"),
    ...
){
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(RColorBrewer)
  
  reduction <- reduction %||% DefaultDimReduc(seurat_object)
  legend.orientation <- match.arg(legend.orientation)
  
  if (merge_legend) {
    if (merge.leg.pos %in% c("right","left")) {
      legend.orientation <- "vertical"
      legend.title.pos   <- "left"
      legend.title.angle <- 90
    } else if (merge.leg.pos %in% c("bottom","top")) {
      legend.orientation <- "horizontal"
      legend.title.pos   <- "top"
      legend.title.angle <- 0
    }
  }
  
  # Handle colors
  if(!is.null(colors_use)){
    if(length(colors_use)<2) stop("colors_use must contain at least 2 colors")
    if(reverse_colors) colors_use <- rev(colors_use)
  } else if (isTRUE(use_viridis)){
    require(viridis)
    colors_use <- viridis::viridis(9, option = viridis_option)
    if(reverse_colors) colors_use <- rev(colors_use)
  } else {
    custom_palettes <- list(
      hotspot = c("navy", "lightblue", "yellow", "orange", "red"),
      five_rainbow = c("blue", "green", "yellow", "orange", "red")
    )
    if(theme_color %in% names(custom_palettes)){
      colors_use <- custom_palettes[[theme_color]]
    } else {
      colors_use <- RColorBrewer::brewer.pal(9, theme_color)
    }
    if(reverse_colors) colors_use <- rev(colors_use)
  }
  
  # Auto pt.size
  if(is.null(pt.size)){
    n_cells <- ncol(seurat_object)
    raster <- raster %||% (n_cells > 2e5)
    pt.size <- if(raster) 1 else min(1583 / n_cells,1)
  }
  
  # Check features
  available_features <- c(rownames(seurat_object), colnames(seurat_object@meta.data))
  all_found_features <- intersect(features, available_features)
  if(length(all_found_features)==0) stop("No features found in Seurat object")
  
  # Blend plot (leave as before but make consistent legend if requested) 
  if(blend){
    if(length(all_found_features)!=2) stop("blend=TRUE requires exactly 2 features")
    p <- FeaturePlot(
      seurat_object,
      features = all_found_features,
      blend = TRUE,
      pt.size = pt.size,
      reduction = reduction,
      raster = raster,
      raster.dpi = raster.dpi,
      label = label,
      ...
    )
    if(no_axes) p <- p & NoAxes()
    if(merge_legend){
      p <- p + plot_layout(guides="collect") & theme(legend.position = "none")
    }
    return(p)
  }
  
  # Compute global scale limits across all features 
  # Use Seurat::FetchData to fetch values (works for both assay features and meta.data)
  # If FetchData fails for features (rare), fallback to GetAssayData.
  fetch_ok <- TRUE
  vals_df <- tryCatch({
    Seurat::FetchData(seurat_object, vars = all_found_features, slot = layer)
  }, error = function(e){
    fetch_ok <<- FALSE
    NULL
  })
  if(!fetch_ok){
    # fallback: try to get from default assay
    assay_name <- Seurat::DefaultAssay(seurat_object)
    vals_df <- do.call(cbind, lapply(all_found_features, function(f){
      mm <- tryCatch(Seurat::GetAssayData(seurat_object[[assay_name]], slot = layer)[f, , drop = TRUE],
                     error = function(e) rep(NA, ncol(seurat_object)))
      return(mm)
    }))
    colnames(vals_df) <- all_found_features
    vals_df <- as.data.frame(vals_df)
  }
  
  # compute global max (ignore NA). Lower bound will be na_cutoff
  global_max <- suppressWarnings(max(as.matrix(vals_df), na.rm = TRUE))
  if(is.infinite(global_max) || is.na(global_max)) global_max <- NA_real_
  color_limits <- c(na_cutoff, global_max)
  # If all values are <= na_cutoff or NA, set an upper slightly > na_cutoff to allow scale to draw
  if(!is.na(global_max) && global_max <= na_cutoff){
    color_limits[2] <- na_cutoff + 1e-6
  }
  
  # Name for the legend (consistent across plots)
  # legend_name <- "Expression level"
  
  # Create single plots (combine = FALSE) and enforce identical scale 
  plot_list <- lapply(all_found_features, function(feat){
    gglist <- FeaturePlot(
      seurat_object,
      features = feat,
      pt.size = pt.size,
      order = order,
      reduction = reduction,
      raster = raster,
      raster.dpi = raster.dpi,
      label = label,
      combine = FALSE,   # crucial: return a list of ggplot objects
      ...
    )
    # Extract ggplot object safely (FeaturePlot returns a list)
    p <- gglist[[1]]
    
    # Apply identical color scale and theme (use + to add ggplot scale)
    p <- p +
      scale_color_gradientn(
        colours = colors_use,
        limits = color_limits,
        na.value = na_color,
        name = NULL,
        guide = guide_colorbar(
          title = NULL,
          title.theme = element_text(angle = legend.title.angle, vjust = 0.5),
          title.position = legend.title.pos,
          direction = legend.orientation,
          barwidth = as.numeric(legend.width),
          barheight = as.numeric(legend.height),
          frame.colour = "black",size=0.2,
          ticks.colour = "black"
        )
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = base_size, face = "bold.italic"),
        legend.title = element_text(size = legend.text.size, face = "bold", hjust = 0.5),
        legend.text = element_text(size = legend.text.size)
      )
    
    if(no_axes) p <- p & NoAxes()
    return(p)
  })
  
  # Combine plots and merge legend 
  plt <- wrap_plots(plot_list, ncol = ncol)
  
  if(merge_legend){
    # Use plot_layout(guides="collect") at patchwork level
    # then set overall legend position
    plt <- plt + plot_layout(guides = "collect") &
      theme(legend.position = merge.leg.pos, legend.justification = "center")
  }
  
  return(plt)
}

