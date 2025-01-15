#SeqArchRplus pipeline
library(seqArchRplus)
library(ChIPseeker)
library(Biostrings)
library(CAGEr)
library(ggplot2)
library(egg)
library(ggpubr)
library(magick)
library(devtools)
library(BSgenome.MorexV3.Gatersleben)
source("seqArchR_utils.R") #slightly edited utility functions from seqArchR+ for our purposes
source("tau_plots_seqArchR.R") #slightly edited seqArchR+ code to work with our tissue specificity data

library(GenomicFeatures)

#setwd("")

# Set seed for reproducibility
set.seed(1234)
genome <- BSgenome.MorexV3.Gatersleben

#for help with
#??seqArchRplus::curate_clusters

#for vignette
#browseVignettes(package = "seqArchRplus")

#get the annotations as a txdb object
txdb <-makeTxDbFromGFF("Hv_Morex.pgsb.Jul2020.gff3", format="gff3")
txdb


#-------------------------------------------------------------------------------------------------------------------------------
## Suppose `merged_sample_names` holds the sample names

merged_sample_names <- c("consensus")
use_dir <- "seqArchR_result"
use_aggl <- "ward.D"
use_dist <- "euclid"
results <- NULL
seqArchR_clusts <- NULL
seqArchRresult <- NULL

for(i in c(1:length(merged_sample_names))){
  sn <- merged_sample_names[i]
  
  #seqArchRresult can be either taken from previous seqArchR run or if one saves the result as RDS file -> loaded
  #seqArchRresult[[i]] <- readRDS(file = paste0("/PATH/TO/SEQARCHR/RESULT/seqArchR_result.rds"))
  results[[sn]] <- seqArchRresult[[i]]
  
  ## get seqArchR clusters custom curated
  seqArchR_clusts[[sn]] <- seqArchRplus::curate_clusters(
    sname = sn,
    use_aggl = use_aggl, use_dist = use_dist,
    seqArchR_result = results[[sn]], iter = 5,
    use_cutk = 2, final = FALSE, dir_path = use_dir
  )
  
}

## Suppose we identify K = 15 upon visual inspection.
#FOR SECONDARY: after inspection this was set to 7 and no manual changes were selected
#FOR PRIMARY: this was set to 15 and the clusters were manually collated as commented bellow
set_cutk <- 15


## Form the lists need_change and change_to for re-assignments
# need_change <- list(c(8, 10, 16))
# change_to <- list(11)

need_change <- list()
change_to <- list()

#9clusters (primary consensus)
#need_change <- list(c(48,49), c(3,5,7), c(14, 25, 8, 9, 11, 13, 1, 2, 10))
#change_to <- list(33, 17, 18)


#only with the "final = TRUE" it shows the changes and final clusters!
## Satisfied> now set final = TRUE
seqArchR_clusts[[sn]] <- seqArchRplus::curate_clusters(sname = sn, 
                                                       use_aggl = use_aggl, use_dist = use_dist,
                                                       seqArchR_result = results[[sn]], iter = 5, 
                                                       
                                                       use_cutk = set_cutk, #***  
                                                       
                                                       need_change = need_change, 
                                                       change_to = change_to, 
                                                       
                                                       final = TRUE, #*** 
                                                       
                                                       dir_path = use_dir)
## This fetches us clusters with custom/curated collation in _arbitrary_ order
## See the next function call to order_clusters_iqw that orders these clusters
## by their median/mean IQW
#---------------------------------------------------------------------------------------------------------


#Generating individual plots via individual functions
#---------------------------------------------------------------------------------------------------------
## Need a large text size because plots are large with high DPI value 
## for production quality
use_txt_size <- 155
merged_sample_names <- c("primary")
proc_data_path <- "/PATH/TO/BEDFILES/"
fname_prefix <- ""
fname_suffix <- "_promoters.bed"
#sample_name <- "consensus" #(previously used for loading BED)
path_to_tau <- "/PATH/TO/TAU/"
tau_names <- NULL
#for primary c(1:9) and paste0("Tau__TC_sample_consensus_cluster", i, ".bed")
#for secondary c(1:7) and paste0("Tau_consensus_alt_cluster", i, ".bed")
for(i in c(1:9)){
  tau_names <- cbind(tau_names, paste0("Tau__TC_sample_consensus_cluster", i, ".bed"))
}

iqw_tpm_pl <- 
  annotations_oneplot_pl <- 
  annotations_list_pl <- 
  seqlogos_oneplot_pl <- 
  seqlogos_list_pl <-  
  stranded_seqlogos_pl <- 
  per_cl_strand_pl <- 
  vector("list", length(merged_sample_names))

names(iqw_tpm_pl) <- 
  names(annotations_oneplot_pl) <- 
  names(annotations_list_pl) <- 
  names(seqlogos_list_pl) <- 
  names(seqlogos_oneplot_pl) <- 
  names(stranded_seqlogos_pl) <- 
  names(per_cl_strand_pl) <- 
  merged_sample_names

seqArchR_clusts_ord <- NULL

for(sn in merged_sample_names){
  
  #sn <- merged_sample_names[3]
  ## Adjust the file path as suitable 
  bed_info_fname <- file.path(proc_data_path, 
                              paste0(fname_prefix, 
                                     sn, fname_suffix))
  
  ##
  info_df <- read.delim(
    file = bed_info_fname,
    sep = "\t", header = TRUE,
    col.names = c(
      "chr", "starts", "ends", "geneId",
      "score", "strand",
      "dominant_ctss", "type", "IQW", "widths", "domTPM", 
      "q_0.1", "q_0.9"
    )
  )
  
  tc_gr <- makeGRangesFromDataFrame(info_df,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames", "seqname",
                                                     "chromosome", "chrom",
                                                     "chr", "chromosome_name",
                                                     "seqid"),
                                    start.field=c("dominant_ctss"),
                                    end.field=c("dominant_ctss"),
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)
  
  
  ## get clusters ordered by median IQW values
  seqArchR_clusts_ord[[sn]] <- seqArchRplus::order_clusters_iqw(
    sname = sn, clusts = seqArchR_clusts[[sn]]$clusters_list,
    info_df = info_df, order_by_median = TRUE)
  ##
  use_clusts <- seqArchR_clusts_ord[[sn]]
  ##
  seqArchRplus::seqs_acgt_image(sname = sn,
                                seqs = results[[sn]]$rawSeqs,
                                seqs_ord = unlist(use_clusts),
                                pos_lab = -50:49, dir_path = use_dir)
  ##
  seqArchRplus::write_seqArchR_cluster_track_bed(sname = sn,
                                                 use_q_bound = F,
                                                 use_as_names = "geneId",
                                                 clusts = use_clusts,
                                                 info_df = info_df,
                                                 one_zip_all = TRUE,
                                                 org_name = fname_prefix,
                                                 dir_path = use_dir,
                                                 include_in_report = FALSE,
                                                 strand_sep = FALSE)
  ##
  iqw_tpm_pl[[sn]] <- iqw_tpm_plots(sname = sn,
                                    dir_path = use_dir,
                                    info_df = info_df,
                                    tau_path = path_to_tau,
                                    tau_names = tau_names,
                                    iqw = TRUE, tpm = TRUE, cons = TRUE, #cons FALSE if no tau plots needed
                                    clusts = use_clusts, 
                                    txt_size = use_txt_size)
  
  
  # iqw_tpm_pl[[sn]]
  # save_plot("IQW_plot_solo.pdf", iqw_tpm_pl[[sn]], ncol = 30, nrow = 10, limitsize = FALSE)
  # 
  
  #Not using seqArchRplus::per_cluster_annotations because we are using function from source(seqArchRplus_utils.R)
  #better for the composed plots (removing some useless labels)
  annotations_oneplot_pl[[sn]] <-
    per_cluster_annotations(
      sname = sn,
      clusts = use_clusts,
      tc_gr = tc_gr,
      cager_obj = NULL,
      qLow = 0.1, qUp = 0.9,
      txdb_obj = txdb,
      tss_region = c(-500,100),
      orgdb_obj = NULL, dir_path = use_dir,
      one_plot = TRUE,
      txt_size = use_txt_size)
  
  # annotations_list_pl[[sn]] <-
  #   seqArchRplus::per_cluster_annotations(
  #     sname = sn,
  #     clusts = use_clusts,
  #     tc_gr = beds[[sample_name]],
  #     cager_obj = NULL,
  #     qLow = 0.1, qUp = 0.9,
  #     txdb_obj = txdb,
  #     tss_region = c(-500,100),
  #     orgdb_obj = NULL, dir_path = use_dir,
  #     one_plot = FALSE,
  #     txt_size = 12)
  
  seqlogos_oneplot_pl[[sn]] <- 
    seqArchRplus::per_cluster_seqlogos(
      sname = sn,
      seqs = results[[sn]]$rawSeqs,
      clusts = use_clusts,
      pos_lab = c(-50:-1, 1:50), bits_yax = "max",
      strand_sep = FALSE, one_plot = TRUE,
      dir_path = use_dir, 
      txt_size = 125)
  
  ## for combining later
  # seqlogos_list_pl[[sn]] <- 
  #   seqArchRplus::per_cluster_seqlogos(sname = sn,
  #                                      seqs = results[[sn]]$rawSeqs,
  #                                      clusts = use_clusts,
  #                                      pos_lab = -50:49, bits_yax = "max",
  #                                      strand_sep = FALSE, one_plot = FALSE,
  #                                      dir_path = use_dir, 
  #                                      txt_size = 12)
  
  stranded_seqlogos_pl[[sn]] <- 
    seqArchRplus::per_cluster_seqlogos(sname = sn,
                                       seqs = results[[sn]]$rawSeqs,
                                       clusts = use_clusts,
                                       pos_lab = -50:49, bits_yax = "max",
                                       info_df = info_df,
                                       strand_sep = TRUE, one_plot = FALSE,
                                       dir_path = use_dir, 
                                       txt_size = 12)
  
  ## Get larger flank FATSA sequences (larger than those used for seqArchR)
  fname <- file.path(proc_data_path,
                     paste0("seqArchR_fasta_", sn,".fa"))
  use_seqs <- Biostrings::readDNAStringSet(filepath = fname, 
                                           format = "FASTA", use.names = TRUE)
  
  plot_motif_heatmaps2(sname = sn, seqs = use_seqs,
                       flanks = c(50),
                       clusts = use_clusts,
                       #BREd, BREu, Y-patch, TATA-box, DPE, MTE (1st ver.), MTE (2nd ver.)
                       #motifs = c("RTDKKKK", "SSRCGCC", "CYTCYYCCYC"),# "RTDKKKK", "SSRCGCC", "CYTCYYCCYC", "TATAWR", "RGWYV", "CSARCSSAACGS", "CGANC"), #motifs = c("GC", "Y", "TATAWAAR", "CAAT", "CSARC", "CGCC", "AACG", "CGTG"), #motifs = c("WW", "SS", "TATAA", "CG", "Y"),motifs = c("TATAAA", "TA", "TATATA", "WW", "SS", "CG", "Y", "GBBRDNHGG"),#motifs = c("RTDKKKK", "Y", "TATAWAAR", "RGWYV", "CSARCSSAACGS"), #motifs = c("WW", "SS", "TATAA", "CG", "Y"), 
                       motifs = c("TA", "CGCC", "WW", "SS"),
                       hm_scale_factor = 0.56, #0.56 - for 4 motifs, 0.75 - for 3 motifs
                       dir_path = use_dir)
  
  pair_colrs <- RColorBrewer::brewer.pal(n = 5, name = "Set3")
  per_cl_strand_pl[[sn]] <- seqArchRplus::per_cluster_strand_dist(sname = sn,
                                                                  clusts = use_clusts,
                                                                  info_df = info_df,
                                                                  dir_path = use_dir,
                                                                  colrs = pair_colrs[4:5])
  
  
}
#------------------------------------------------------------------------------------------------------------------------

#GO MANUALLY THROUGH ALL SAMPLES, SET HEIGHTS ETC., FOR THE BEST VISUALS!!!!!!!
library(cowplot)
sn <- merged_sample_names[1]

p <- grid.arrange(arrangeGrob(iqw_tpm_pl[[sn]], grid::nullGrob(),
                              heights = c(1, 0)),
                  arrangeGrob(grid::nullGrob(), seqlogos_oneplot_pl[[sn]], grid::nullGrob(),
                              heights = c(0.03, 0.85, 0.12)),
                  annotations_oneplot_pl[[sn]],
                  ncol = 3)
save_plot("/PATH/TO/RESULTS/seqArchR_primary_plot.pdf", p, ncol = 30, nrow = 17, limitsize = FALSE)
