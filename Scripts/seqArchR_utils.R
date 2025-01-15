## Functions for downstream analysis of seqArchR clusters
##
## Each of these functions performs the task for one specified sample
## 1. Tag clusters as bed files for viewing in Genome browser/IGV DONE
## 2. Per cluster annotation percentages DONE
## 3. Arranging clusters by IQW, TPM for comparison of architectures DONE
## 4. Seqlogos, usual and strand-separated.
## 5. How to visualize prevalence of clusters on different chromosomes
## 6. Heatmaps arranged by IQW
## 7. Manual curation of clusters

## Handle multicore bpparam object creation using BiocParallel
##
handle_multicore <- function(crs = 1, parallelize = FALSE) {

  #### Start cluster only once -- using BiocParallel
  if (parallelize) {
    if (.Platform$OS.type == "windows") {
      if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
      cl <- BiocParallel::SnowParam(
        workers = crs, tasks = crs,
        exportglobals = TRUE
      )
      # cli::cli_alert_info("Parallelization: Yes")
    } else {
      if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
      cl <- BiocParallel::MulticoreParam(workers = crs, tasks = crs)
      # cli::cli_alert_info("Parallelization: Yes")
    }
  } else {
    cl <- BiocParallel::SerialParam()
    # cli::cli_alert_info("Parallelization: No")
  }
  return(cl)
}
## =============================================================================



get_strand_plot_title <- function(this_id, nclust, this_n = NULL,
                                   tot_n, strand_val = "+") {
  if (!is.null(strand_val)) {
    if (is.null(this_n)) {
      stop("'this_n' cannot be NULL when strand_val is not NULL")
    }
  }

  clust_names <- sort(as.character(seq(nclust)))

  if (is.null(strand_val)) {
    title_str <- paste0(
      "(", this_id, "/", nclust,
      ") Arch `", clust_names[this_id], "': ",
      tot_n, " sequences"
    )
  } else {
    if (strand_val == "+") {
      title_str <- paste0(
        "(", this_id, "/", nclust,
        ") Arch `", clust_names[this_id], "': ",
        this_n, " (+ strand) /", tot_n
      )
    } else if (strand_val == "-") {
      title_str <- paste0(
        "(", this_id, "/", nclust,
        ") Arch `", clust_names[this_id], "': ",
        this_n, " (- strand) /", tot_n
      )
    }
  }


  title_str
}
## =============================================================================

get_strand_specific_indices <- function(df, seq_ids_in_cl, strand_val = "+") {
  return(seq_ids_in_cl[which(df$strand[seq_ids_in_cl] == strand_val)])
}
## =============================================================================

## Handle per sample result directory
handle_per_sample_result_dir <- function(sname, dir_path) {
  result_dir_path <- file.path(dir_path, paste0(sname, "_results"))
  stopifnot(check_and_create_dir(result_dir_path))
  result_dir_path
}
## =============================================================================

## check_and_create_dir
check_and_create_dir <- function(dir_path) {
  creation_ok <- FALSE
  if (!dir.exists(dir_path)) {
    cli::cli_alert_warning(paste0("Creating directory: ", dir_path))
    creation_ok <- dir.create(dir_path)
  } else {
    cli::cli_alert_warning(paste0("Directory exists: ", dir_path))
    creation_ok <- TRUE
  }
  ## TRUE: success in creation; FALSE: otherwise
  return(creation_ok)
}
## =============================================================================





write_empty_string <- function() {
  cat_str <- paste0("<a href= >Empty", "</a>\n")
  cat_str
}
## =============================================================================



get_clust_id_column <- function(info_df, clusts, use_prefix, use_suffix) {
  ## Add new column noting the cluster IDs from seqArchR result
  clust_lab <- rep("0", length(unlist(clusts)))
  ##
  clust_labels <- make_cluster_labels(clusts, use_prefix, use_suffix)

  for (i in seq_along(clusts)) {
    clust_lab[clusts[[i]]] <- clust_labels[i]
  }

  ## Return as a factor because this is directly used as y-axis labels
  ## Reverse the levels order so that the clusters are placed top-to-botom
  clust_lab <- factor(clust_lab, levels = rev(clust_labels))

  clust_lab
}
## =============================================================================

make_cluster_labels <- function(clust, use_prefix, use_suffix) {
  clust_lens <- lengths(clust)
  clust_labels <- paste(paste0("(", clust_lens, ") "),
                        use_prefix, seq_along(clust), use_suffix,
                        sep = ""
  )
  clust_labels
}
## =============================================================================


## Use this function to get n colors (in sequence) from the specified palette
## -- Can specify n greater than that existing in a palette, in which case
## the colors are recycled and returned
## -- Can also specify a list of colors of choice < n, which are recycled and
## returned. Useful when random colors from a palette are required
get_ncolors <- function(n, palname = "Set1", clrs = NULL) {

  ## Recycle color vector from RColorBrewer
  use_colors <- clrs
  if (is.null(clrs)) {
    use_colors <- RColorBrewer::brewer.pal(n = n, name = palname)
  }

  nColor <- length(use_colors)
  if (n <= nColor) {
    n_colors <- use_colors[seq_len(n)]
    return(n_colors)
  }

  rep_times <- base::ceiling((n - nColor) / nColor)
  if (n %% nColor == 0) rep_times <- rep_times + 1
  additional <- ((n - nColor) %% nColor)
  col_idx <- c(rep(seq_len(nColor), rep_times), seq_len(additional))
  n_colors <- use_colors[col_idx]
  n_colors
}
## =============================================================================


## Handles preparing a tagcluster GRanges object from either a CAGEexp object
## if that is provided.
## Also, it can be prepared using a BED file provided.
handle_tc_cager <- function(tc_gr = NULL, cager_obj = NULL, sname,
                             qLow = 0.1, qUp = 0.9) {
  if (is.null(tc_gr)) {
    if (is.null(cager_obj)) {
      stop(
        "`tc_gr` is NULL, did you forgot to supply the `cager_obj`?",
        "Please also specify `qLow` and `qUp` with the `cager_obj`"
      )
    } else {
      if (!requireNamespace("CAGEr", quietly = TRUE)) {
        stop(
          "Please install R package 'CAGEr' to automatically ",
          "extract CAGE tag clusters."
        )
      }
      tc_gr <- get_cage_tc_cager(
        cager_obj = cager_obj, sname = sname,
        qLow = qLow, qUp = qUp
      )
      seqlevels(tc_gr) <- GenomeInfoDb::seqlevels(
        CAGEr::CTSStagCountSE(cager_obj)
      )
      seqinfo(tc_gr) <- GenomeInfoDb::seqinfo(
        CAGEr::CTSStagCountSE(cager_obj)
      )
    }
  }

  if(!(is(tc_gr, "GRanges") || is(tc_gr, "TagClusters"))){
    if(.check_bedfile(tc_gr)){
      cli::cli_alert_info("Reading from Bed file")
      df <- read.delim(tc_gr, header = TRUE, sep="\t")
      tc_gr <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                       keep.extra.columns = TRUE)
    }
    return(list(tc_gr, "bed"))
  }else{
    return(list(tc_gr, NULL))
  }
  stop("Expecting `tc_gr` to be a GRanges object or a bedfile")

}
## =============================================================================

check_bedfile <- function(bedfname){
  ## Just check if bedfname is valid and exists
  if(file.exists(bedfname))
    return(TRUE)
  return(FALSE)
}
## =============================================================================

get_cage_tc_cager <- function(cager_obj, sname, qLow, qUp) {
  any_null <- any(unlist(lapply(list(qLow, qUp), is.null)))
  if (any_null) stop("Please specify both `qLow` and `qUp`.")
  message("Using qLow = ", qLow, " and qUp = ", qUp)
  tc_gr <- CAGEr::tagClustersGR(cager_obj,
                                sample = sname,
                                returnInterquantileWidth = TRUE,
                                qLow = qLow, qUp = qUp
  )
  tc_gr
}
## =============================================================================
per_cluster_annotations <- function(sname = NULL, clusts = NULL,
                                    tc_gr = NULL,
                                    cager_obj = NULL,
                                    qLow = 0.1,
                                    qUp = 0.9,
                                    txdb_obj = NULL,
                                    tss_region = NULL,
                                    orgdb_obj = NULL,
                                    one_plot = TRUE,
                                    dir_path = NULL,
                                    txt_size = 12,
                                    use_suffix = NULL, use_prefix = "C",
                                    n_cores = 1) {
  
  cli::cli_h1(paste0("All clusters' genomic annotations"))
  cli::cli_h2(paste0("Sample: ", sname))
  ## Check all needed arguments supplied
  
  
  ## Prepare tc_gr
  tc_gr2 <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
  stopifnot(!is.null(tc_gr2))
  ## If tc_gr was prepared from a BED file, populate the clusts arg
  
  if(!is.null(tc_gr2[[2]]) && tc_gr2[[2]] == "bed"){
    clusts <- seq(length(tc_gr2[[1]]))
  }
  
  ## clusts should be a list
  if(!is.list(clusts)) clusts <- list(clusts)
  ##
  tc_gr <- tc_gr2[[1]]
  
  ##
  ## as many records in tc_gr as number of sequence IDs in clusts
  # if(length(tc_gr) != sum(lengths(clusts))){
  #     stop("Nb. of records in `tc_gr` should match the nb. of sequence IDs
  #         in `clusts`")
  # }
  ##
  parallelize <- FALSE
  if (n_cores > 1) parallelize <- TRUE
  if(parallelize)
    bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
  ##
  if(!is.null(dir_path)){
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    fname <- file.path(result_dir_path,
                       paste0("Clusterwise_annotations.pdf"))
  }
  ##
  clust_labels <- .make_cluster_labels(clust = clusts, use_prefix, use_suffix)
  
  ## Without BiocParallel
  # clustwise_anno <- lapply(clusts, function(x) {
  #     foo_anno <- ChIPseeker::annotatePeak(tc_gr[x, ],
  #         tssRegion = tss_region,
  #         TxDb = txdb_obj,
  #         annoDb = orgdb_obj
  #     )
  #     foo_anno
  # })
  # names(clustwise_anno) <- seq(1, length(clusts))
  # return(clustwise_anno) ## to check
  # ## without dplyr?
  # df_list <- lapply(clustwise_anno, function(x) x@annoStat)
  
  ##
  peakAnno <- ChIPseeker::annotatePeak(tc_gr,
                                       tssRegion = tss_region,
                                       TxDb = txdb_obj,
                                       annoDb = orgdb_obj
  )
  ##
  
  clustwise_anno <- lapply(clusts, function(x){
    .get_anno_stat(peakAnno, x)
  })
  names(clustwise_anno) <- seq(1, length(clusts))
  ##
  
  colrs <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
  ## Note: Generally, speed is not an issue when the number of clusters is
  ## around 40-50 or even up to 100.
  ## So we can use the base::rbind instead of the dplyr::bind_rows which is
  ## generally faster.
  ##
  
  sam <- do.call("rbind", clustwise_anno)
  per_df_nrows <- unlist(lapply(clustwise_anno, nrow))
  ## Add a column clust that we need downstream
  # sam$clust <- rep(seq_along(clustwise_anno), times = per_df_nrows)
  
  sam$clust <- rep(clust_labels, times = per_df_nrows)
  sam$clust <- factor(sam$clust, levels = rev(clust_labels))
  
  sam <- sam[, c(3, 1, 2)]
  rownames(sam) <- seq_len(sum(per_df_nrows))
  ##
  
  sam$Feature <- factor(sam$Feature,
                        levels = c(
                          "Promoter", "5' UTR", "1st Exon", "Other Exon",
                          "1st Intron", "Other Intron", "Downstream (<=300)",
                          "3' UTR", "Distal Intergenic"
                        )
  )
  ##
  names(colrs) <- levels(sam$Feature)
  ##
  if (one_plot) {
    ## Solution using custom plotting
    clustwise_annobar <- .get_prop_anno_oneplot(sam,
                                                txt_size = txt_size,
                                                colrs = colrs
    )
    
    ## Put in the cluster names when printing independently
    ind_print_plot <- clustwise_annobar +
      ggplot2::scale_y_discrete(
        labels = rev(clust_labels),
        expand = expansion(add = c(0.5, 0.5))
      ) +
      ggplot2::theme(axis.text.y = element_blank())
    ind_print_plot <- ind_print_plot +
      ggplot2::labs(x = "Proportion", y = "")
    ##
    if(!is.null(dir_path)){
      ggplot2::ggsave(
        filename = fname, plot = ind_print_plot,
        device = "pdf", width = 15, height = 20, units = "in",
        dpi = 300
      )
    }
    return(ind_print_plot)
  } else {
    ## return individual plots as a list
    sam_split <- split(sam, f = sam$clust)
    
    ## Without BiocParallel
    annobar_list <- lapply(seq_along(sam_split), function(x) {
      ##
      anno_df <- sam_split[[x]]
      use_clust_label <- clust_labels[x]
      anno_df$Feature <- factor(anno_df$Feature,
                                levels = levels(sam$Feature)
      )
      anno_pl <- .get_prop_anno_listplot(anno_df,
                                         txt_size = txt_size,
                                         colrs = colrs
      )
      anno_pl <- anno_pl +
        ggplot2::labs(x = "Proportion", y = "Clusters")
      
    })
    return(annobar_list)
  }
}
## =============================================================================


## anno_df holds information on the annoStats from the ChIPSeeker annotatePeak
## return object
## Expected columns are 'Frequency', 'Feature', 'clust'.
##
##
.get_prop_anno_oneplot <- function(anno_df, txt_size, colrs) {
  clustwise_annobar <- ggplot(
    anno_df,
    aes(y = clust, x = Frequency, fill = Feature)
  ) +
    geom_bar(
      stat = "identity", width = 0.6,
      position = position_fill(reverse = TRUE)
    ) +
    ggplot2::scale_fill_manual(name = "", values = colrs) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = element_blank(),
      axis.ticks.length.y.left = unit(0.1, units = "cm"),
      axis.ticks.length.x.bottom = unit(0.1, units = "cm"),
      axis.text.x.bottom = element_text(size = txt_size),
      axis.title.x.bottom = element_text(size = txt_size),
      axis.title.y.left = element_text(size = txt_size),
      legend.position = "bottom",
      legend.text = element_text(size = txt_size)
    ) +
    ggplot2::scale_y_discrete(expand = expansion(add = c(0.75, 0))) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.01)) +
    ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE)) +
    NULL
  ##
  clustwise_annobar
}
## =============================================================================

.get_anno_stat <- function(peakAnno, idx){
  
  anno_terms1 <- as.data.frame(peakAnno)[idx,]$annotation
  anno_terms_splits <- strsplit(anno_terms1, " ")
  anno_terms <- unlist(lapply(anno_terms_splits, function(x){
    if(is.na(x[1])){ "NA";
    }else if(grepl("Intron", x[1]) || grepl("Exon", x[1])){
      if(x[4] == 1){
        paste("1st", x[1])
      }else{
        paste("Other", x[1])
      }
    }else if(grepl("Promoter", x[1])){
      x[1]
    }else{ paste(x[1], x[2]) }
  }))
  ##
  wordFreqTable <- table(anno_terms)
  ##
  wordFreqDF <- data.frame(Feature = rownames(wordFreqTable),
                           Frequency = as.vector(100*wordFreqTable/sum(wordFreqTable)))
  wordFreqDF
}
## =============================================================================

.get_prop_anno_listplot <- function(anno_df, txt_size, colrs) {
  ##
  pl <- ggplot(anno_df, aes(y = clust, x = Frequency, fill = Feature)) +
    # ggplot2::theme_bw() +
    ggplot2::geom_bar(
      stat = "identity", width = 0.9,
      position = position_fill(reverse = TRUE)
    ) +
    ggplot2::scale_fill_manual(values = colrs) +
    ggplot2::scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +
    ggplot2::scale_x_continuous(expand = expansion(mult = 0.05)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y = element_blank(),
      axis.title = element_text(size = txt_size),
      axis.text = element_text(size = txt_size),
      axis.text.y = element_blank(),
      axis.ticks.length.y.left = unit(0.2, units = "cm"),
      panel.grid = element_blank(),
      legend.text = element_text(size = txt_size),
      legend.title = element_text(size = txt_size),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE)) +
    # # ggplot2::labs(x = "Percentage (%)", y = "Cluster") +
    NULL
  return(pl)
}
## =============================================================================


.get_fixed_anno_ord <- function() {
  anno_terms_ord <- c(
    "Promoter", "5' UTR", "3' UTR", "1st Exon",
    "Other Exon", "1st Intron", "Intron",
    "Downstream",
    "Distal Intergenic"
  )
  
  anno_terms_ord
}
## =============================================================================


.get_named_colors <- function(anno_terms_ord, palname = "Set1") {
  use_colors <- .get_ncolors(n = length(anno_terms_ord), palname = palname)
  if (palname == "Set1") use_colors <- use_colors[c(2:length(use_colors), 1)]
  names(use_colors) <- anno_terms_ord
  use_colors
}


## Handle multicore bpparam object creation using BiocParallel
##
.handle_multicore <- function(crs = 1, parallelize = FALSE) {
  
  #### Start cluster only once -- using BiocParallel
  if (parallelize) {
    if (.Platform$OS.type == "windows") {
      if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
      cl <- BiocParallel::SnowParam(
        workers = crs, tasks = crs,
        exportglobals = TRUE
      )
      # cli::cli_alert_info("Parallelization: Yes")
    } else {
      if (is.null(crs)) crs <- BiocParallel::multicoreWorkers()
      cl <- BiocParallel::MulticoreParam(workers = crs, tasks = crs)
      # cli::cli_alert_info("Parallelization: Yes")
    }
  } else {
    cl <- BiocParallel::SerialParam()
    # cli::cli_alert_info("Parallelization: No")
  }
  return(cl)
}
## =============================================================================



.get_strand_plot_title <- function(this_id, nclust, this_n = NULL,
                                   tot_n, strand_val = "+") {
  if (!is.null(strand_val)) {
    if (is.null(this_n)) {
      stop("'this_n' cannot be NULL when strand_val is not NULL")
    }
  }
  
  clust_names <- sort(as.character(seq(nclust)))
  
  if (is.null(strand_val)) {
    title_str <- paste0(
      "(", this_id, "/", nclust,
      ") Arch `", clust_names[this_id], "': ",
      tot_n, " sequences"
    )
  } else {
    if (strand_val == "+") {
      title_str <- paste0(
        "(", this_id, "/", nclust,
        ") Arch `", clust_names[this_id], "': ",
        this_n, " (+ strand) /", tot_n
      )
    } else if (strand_val == "-") {
      title_str <- paste0(
        "(", this_id, "/", nclust,
        ") Arch `", clust_names[this_id], "': ",
        this_n, " (- strand) /", tot_n
      )
    }
  }
  
  
  title_str
}
## =============================================================================

.get_strand_specific_indices <- function(df, seq_ids_in_cl, strand_val = "+") {
  return(seq_ids_in_cl[which(df$strand[seq_ids_in_cl] == strand_val)])
}
## =============================================================================

## Handle per sample result directory
.handle_per_sample_result_dir <- function(sname, dir_path) {
  result_dir_path <- file.path(dir_path, paste0(sname, "_results"))
  stopifnot(.check_and_create_dir(result_dir_path))
  result_dir_path
}
## =============================================================================

## check_and_create_dir
.check_and_create_dir <- function(dir_path) {
  creation_ok <- FALSE
  if (!dir.exists(dir_path)) {
    cli::cli_alert_warning(paste0("Creating directory: ", dir_path))
    creation_ok <- dir.create(dir_path)
  } else {
    cli::cli_alert_warning(paste0("Directory exists: ", dir_path))
    creation_ok <- TRUE
  }
  ## TRUE: success in creation; FALSE: otherwise
  return(creation_ok)
}
## =============================================================================





.write_empty_string <- function() {
  cat_str <- paste0("<a href= >Empty", "</a>\n")
  cat_str
}
## =============================================================================



.get_clust_id_column <- function(info_df, clusts, use_prefix, use_suffix) {
  ## Add new column noting the cluster IDs from seqArchR result
  clust_lab <- rep("0", length(unlist(clusts)))
  ##
  clust_labels <- .make_cluster_labels(clusts, use_prefix, use_suffix)
  
  for (i in seq_along(clusts)) {
    clust_lab[clusts[[i]]] <- clust_labels[i]
  }
  
  ## Return as a factor because this is directly used as y-axis labels
  ## Reverse the levels order so that the clusters are placed top-to-botom
  clust_lab <- factor(clust_lab, levels = rev(clust_labels))
  
  clust_lab
}
## =============================================================================

.make_cluster_labels <- function(clust, use_prefix, use_suffix) {
  clust_lens <- lengths(clust)
  clust_labels <- paste(paste0("(", clust_lens, ") "),
                        use_prefix, seq_along(clust), use_suffix,
                        sep = ""
  )
  clust_labels
}
## =============================================================================


## Use this function to get n colors (in sequence) from the specified palette
## -- Can specify n greater than that existing in a palette, in which case
## the colors are recycled and returned
## -- Can also specify a list of colors of choice < n, which are recycled and
## returned. Useful when random colors from a palette are required
.get_ncolors <- function(n, palname = "Set1", clrs = NULL) {
  
  ## Recycle color vector from RColorBrewer
  use_colors <- clrs
  if (is.null(clrs)) {
    use_colors <- RColorBrewer::brewer.pal(n = n, name = palname)
  }
  
  nColor <- length(use_colors)
  if (n <= nColor) {
    n_colors <- use_colors[seq_len(n)]
    return(n_colors)
  }
  
  rep_times <- base::ceiling((n - nColor) / nColor)
  if (n %% nColor == 0) rep_times <- rep_times + 1
  additional <- ((n - nColor) %% nColor)
  col_idx <- c(rep(seq_len(nColor), rep_times), seq_len(additional))
  n_colors <- use_colors[col_idx]
  n_colors
}
## =============================================================================


## Handles preparing a tagcluster GRanges object from either a CAGEexp object
## if that is provided.
## Also, it can be prepared using a BED file provided.
.handle_tc_cager <- function(tc_gr = NULL, cager_obj = NULL, sname,
                             qLow = 0.1, qUp = 0.9) {
  if (is.null(tc_gr)) {
    if (is.null(cager_obj)) {
      stop(
        "`tc_gr` is NULL, did you forgot to supply the `cager_obj`?",
        "Please also specify `qLow` and `qUp` with the `cager_obj`"
      )
    } else {
      if (!requireNamespace("CAGEr", quietly = TRUE)) {
        stop(
          "Please install R package 'CAGEr' to automatically ",
          "extract CAGE tag clusters."
        )
      }
      tc_gr <- .get_cage_tc_cager(
        cager_obj = cager_obj, sname = sname,
        qLow = qLow, qUp = qUp
      )
      seqlevels(tc_gr) <- GenomeInfoDb::seqlevels(
        CAGEr::CTSStagCountSE(cager_obj)
      )
      seqinfo(tc_gr) <- GenomeInfoDb::seqinfo(
        CAGEr::CTSStagCountSE(cager_obj)
      )
    }
  }
  
  if(!(is(tc_gr, "GRanges") || is(tc_gr, "TagClusters"))){
    if(.check_bedfile(tc_gr)){
      cli::cli_alert_info("Reading from Bed file")
      df <- read.delim(tc_gr, header = TRUE, sep="\t")
      tc_gr <- GenomicRanges::makeGRangesFromDataFrame(df,
                                                       keep.extra.columns = TRUE)
    }
    return(list(tc_gr, "bed"))
  }else{
    return(list(tc_gr, NULL))
  }
  stop("Expecting `tc_gr` to be a GRanges object or a bedfile")
  
}
## =============================================================================

.check_bedfile <- function(bedfname){
  ## Just check if bedfname is valid and exists
  if(file.exists(bedfname))
    return(TRUE)
  return(FALSE)
}
## =============================================================================

.get_cage_tc_cager <- function(cager_obj, sname, qLow, qUp) {
  any_null <- any(unlist(lapply(list(qLow, qUp), is.null)))
  if (any_null) stop("Please specify both `qLow` and `qUp`.")
  message("Using qLow = ", qLow, " and qUp = ", qUp)
  tc_gr <- CAGEr::tagClustersGR(cager_obj,
                                sample = sname,
                                returnInterquantileWidth = TRUE,
                                qLow = qLow, qUp = qUp
  )
  tc_gr
}
## =============================================================================



