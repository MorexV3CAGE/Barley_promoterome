#SeqArchR pipeline
library(seqArchR)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.MorexV3.Gatersleben)
library(dplyr)

# Set seed for reproducibility
set.seed(1234)

#set BSgenome into variable
genome <- BSgenome.MorexV3.Gatersleben

#---------------------------------------------------------------------------
#function to run seqArchR for any input and have an output of RDS files
run_seqArchR <- function(names, path_to_files, strictsness, fifth, output_path){
  genome <- BSgenome.MorexV3.Gatersleben
  promoter_set <- NULL
  granges_pooled <- NULL
  promoter_seqs_pooled <- NULL
  inputSeqs_direct_pooled <- NULL
  nSeqs_pooled <- NULL
  positions_pooled <- NULL
  promoter_seqs_pooled_fasta <- NULL
  
  for(i in names){
    #get input into DNAStringSet
    #---------------------------------------------------------------------------
    #1. Input bed files as dataframes
    #4DAG/8DAP/24DAP sample files
    j <- which(names == i)
    promoter_set[[i]] <- read.table(path_to_files[j], header = TRUE)
    
    #changes so the granges convert doesn't have a problem
    names(promoter_set[[i]])[names(promoter_set[[i]]) == 'start'] <- 'starts'
    names(promoter_set[[i]])[names(promoter_set[[i]]) == 'end'] <- 'ends'
    names(promoter_set[[i]])[names(promoter_set[[i]]) == 'width'] <- 'length'
    
    #filter out ribosomal sequences (should be already done beforehand, so shouldn't change the dataset)
    #promoter_set[[i]] <- subset(promoter_set[[i]], seqnames != "chr5H" & seqnames != "chr6H" | seqnames == "chr5H" & dominantTSS<52608306 | seqnames == "chr5H" & dominantTSS>53499223 | seqnames == "chr6H" & dominantTSS<81918150 | seqnames == "chr6H" & dominantTSS>82454047)
    
    #2. These dataframes are converted to GRanges, with start and end set to the dominant peak, so the extraction can do its thing
    #RETURN TO GRANGES
    granges_pooled[[i]] <- makeGRangesFromDataFrame(promoter_set[[i]],
                                                    keep.extra.columns=TRUE,
                                                    ignore.strand=FALSE,
                                                    seqinfo=NULL,
                                                    seqnames.field=c("seqnames", "seqname",
                                                                     "chromosome", "chrom",
                                                                     "chr", "chromosome_name",
                                                                     "seqid"),
                                                    start.field=c("dominantTSS"),
                                                    end.field=c("dominantTSS"),
                                                    strand.field="strand",
                                                    starts.in.df.are.0based=FALSE)
    
    #consensus_output <- as.data.frame(granges_pooled)
    #get fasta files (later needed for seqArchR+)
    promoter_seqs_pooled_fasta[[i]] <- granges_pooled[[i]] %>%
      promoters(upstream = 500, downstream = 500) %>%
      getSeq(genome, .)
    
    writeXStringSet(promoter_seqs_pooled_fasta[[i]], paste0(output_path, "seqArchR_fasta_", i, ".fa"))
    
    #3. Get the sequences
    promoter_seqs_pooled[[i]] <- granges_pooled[[i]] %>%
      promoters(upstream = 50, downstream = 50) %>%
      getSeq(genome, .)
    
    
    # Creation of one-hot encoded data matrix from a DNAStringSet object
    inputSeqs_direct_pooled[[i]] <- seqArchR::get_one_hot_encoded_seqs(seqs = promoter_seqs_pooled[[i]], 
                                                                       sinuc_or_dinuc = "dinuc")
    
    nSeqs_pooled[[i]] <- length(promoter_seqs_pooled[[i]])
    
    positions_pooled[[i]] <- seq(1, Biostrings::width(promoter_seqs_pooled[[i]][1]))
    
    #Visualize
    seqArchR::viz_seqs_acgt_mat(as.character(promoter_seqs_pooled[[i]]), 
                                pos_lab = positions_pooled[[i]], save_fname = NULL)
    dev.copy(jpeg,paste0(output_path, i, "_heatmap.jpg"),height=2000,width=1600,res=200) 
    dev.off()
    
    # Set seqArchR configuration
    seqArchRconfig <- seqArchR::set_config(
      parallelize = TRUE,
      n_cores = 2,
      n_runs = 100,
      k_min = 1,
      k_max = 20,
      mod_sel_type = "stability",
      bound = strictsness, #for MorexV3 10^-8 was ideal
      chunk_size = fifth[j], #depends on the number of CAGEr candidates (split by 5 was good)
      result_aggl = "ward.D", 
      result_dist = "euclid",
      flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
                   plot = FALSE)
    )
    
    # Call/Run seqArchR
    seqArchRresult[[i]] <- seqArchR::seqArchR(config = seqArchRconfig,
                                              seqs_ohe_mat = inputSeqs_direct_pooled[[i]],
                                              seqs_raw = promoter_seqs_pooled[[i]],
                                              seqs_pos = positions_pooled[[i]],
                                              total_itr = 5,
                                              set_ocollation = c(FALSE, TRUE, TRUE, TRUE, FALSE))
    
    saveRDS(seqArchRresult[[i]], paste0(output_path, "seqArchR_result_", i, ".rds"))
  }
}


#SEQARCHR IN A SINGLE FUNCTION
#---------------------------------------------------------------------------
#1. Input names you want to save the results under and PATH to the input files
names <- c("primary", "secondary")
path_to_files <- c("/PATH/TO/BEDFILES/primary_seqarchr.bed",
                   "/PATH/TO/BEDFILES/secondaty_seqarchr.bed") #supose that all of the inputs have the same structure as the example file Data_S1.txt
strictsness <- 10^-8 #10^-8 ideal for MorexV3, 10^-6 = less strict, 10^-10 = more strict
fifth <- c(4500, 500) #set as 1/5 of your dataset/s

output_path <- "/PATH/TO/RESULTS/"

run_seqArchR(names, path_to_files, strictsness, fifth, output_path)


#some additional visualizations according to the manual
seqArchR::viz_bas_vec(feat_mat = get_clBasVec_m(seqArchRresult, 1), 
                      ptype = c("heatmap", "seqlogo"), method = "bits", 
                      sinuc_or_dinuc = "dinuc")

seqArchR::viz_bas_vec(feat_mat = get_clBasVec_m(seqArchRresult, 2), 
                      ptype = c("heatmap", "seqlogo"), method = "bits", 
                      sinuc_or_dinuc = "dinuc")




seqArchR::viz_seqs_acgt_mat(seqs_str(seqArchRresult, iter = 2, ord = TRUE),
                            pos_lab = positions)



