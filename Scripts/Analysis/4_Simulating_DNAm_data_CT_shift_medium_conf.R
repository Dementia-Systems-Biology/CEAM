source("Scripts/Functions.R")

# load required packages
check_and_load_libraries(c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "patchwork",
  "purrr",
  "MCMCpack",
  "data.table",
  "parallel",
  "pbapply",
  "msm",
  "truncnorm"
))

load("Lvl2_CpG_sets_Annotated.RData")

# find set sizes of real data
set_sizes <- c(length(NeuN_lvl2_set$cpg), 
               length(IRF8_lvl2_set$cpg), 
               length(SOX10_lvl2_set$cpg), 
               length(TN_lvl2_set$cpg)) # set sizes of lvl2 confidence
names(set_sizes) <- c("NeuN", "IRF8", "SOX10", "TN")
background_size <- 785614

# initialize the nr of CpGs and samples that should be simulated
sample_num <- 250
cell_types_num <- 4
scaling_factor <- 0.25

# Calculate pairwise intersections
NeuN_SOX10_intersect <- intersect(NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg)
NeuN_IRF8_intersect <- intersect(NeuN_lvl2_set$cpg, IRF8_lvl2_set$cpg)
NeuN_TN_intersect <- intersect(NeuN_lvl2_set$cpg, TN_lvl2_set$cpg)
SOX10_IRF8_intersect <- intersect(SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg)
SOX10_TN_intersect <- intersect(SOX10_lvl2_set$cpg, TN_lvl2_set$cpg)
IRF8_TN_intersect <- intersect(IRF8_lvl2_set$cpg, TN_lvl2_set$cpg)

# Calculate three-way intersections using Reduce
NeuN_SOX10_IRF8_intersect <- Reduce(intersect, list(NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg))
NeuN_SOX10_TN_intersect <- Reduce(intersect, list(NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg, TN_lvl2_set$cpg))
NeuN_IRF8_TN_intersect <- Reduce(intersect, list(NeuN_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg))
SOX10_IRF8_TN_intersect <- Reduce(intersect, list(SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg))


sim_set_sizes <- c(
  # Set Sizes (after setdiff operation and scaling)
  NeuN_set_size = round(length(Reduce(setdiff, list(
    NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor),
  
  IRF8_set_size = round(length(Reduce(setdiff, list(
    IRF8_lvl2_set$cpg, SOX10_lvl2_set$cpg, NeuN_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor),
  
  SOX10_set_size = round(length(Reduce(setdiff, list(
    SOX10_lvl2_set$cpg, NeuN_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor),
  
  TN_set_size = round(length(Reduce(setdiff, list(
    TN_lvl2_set$cpg, SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg, NeuN_lvl2_set$cpg
  ))) * scaling_factor),
  
  # Calculate the sizes, subtracting the three-way intersections from the two-way ones
  NeuN_SOX10_size = round((length(NeuN_SOX10_intersect) - length(NeuN_SOX10_IRF8_intersect) - length(NeuN_SOX10_TN_intersect)) * scaling_factor),
  NeuN_IRF8_size = round((length(NeuN_IRF8_intersect) - length(NeuN_SOX10_IRF8_intersect) - length(NeuN_IRF8_TN_intersect)) * scaling_factor),
  NeuN_TN_size = round((length(NeuN_TN_intersect) - length(NeuN_SOX10_TN_intersect) - length(NeuN_IRF8_TN_intersect)) * scaling_factor),
  SOX10_IRF8_size = round((length(SOX10_IRF8_intersect) - length(NeuN_SOX10_IRF8_intersect) - length(SOX10_IRF8_TN_intersect)) * scaling_factor),
  SOX10_TN_size = round((length(SOX10_TN_intersect) - length(NeuN_SOX10_TN_intersect) - length(SOX10_IRF8_TN_intersect)) * scaling_factor),
  IRF8_TN_size = round((length(IRF8_TN_intersect) - length(NeuN_IRF8_TN_intersect) - length(SOX10_IRF8_TN_intersect)) * scaling_factor),
  
  # Three-way Intersections (after intersect operation and scaling)
  NeuN_SOX10_IRF8 = round(length(Reduce(intersect, list(
    NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg
  ))) * scaling_factor),
  
  NeuN_SOX10_TN = round(length(Reduce(intersect, list(
    NeuN_lvl2_set$cpg, SOX10_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor),
  
  NeuN_IRF8_TN = round(length(Reduce(intersect, list(
    NeuN_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor),
  
  SOX10_IRF8_TN = round(length(Reduce(intersect, list(
    SOX10_lvl2_set$cpg, IRF8_lvl2_set$cpg, TN_lvl2_set$cpg
  ))) * scaling_factor)
)

# find how many CpGs are in a set and how many are not
cpg_in_set_num <- as.numeric(cumsum(sim_set_sizes)[length(cumsum(sim_set_sizes))])
cpg_out_set_num <- round((background_size*scaling_factor) - cpg_in_set_num)

# count the total CpGs in each cell type set
NeuN_set_size <- sum(as.numeric(sim_set_sizes[grepl("NeuN",names(sim_set_sizes))]))
IRF8_set_size <- sum(as.numeric(sim_set_sizes[grepl("IRF8",names(sim_set_sizes))]))
SOX10_set_size <- sum(as.numeric(sim_set_sizes[grepl("SOX10",names(sim_set_sizes))]))
TN_set_size <- sum(as.numeric(sim_set_sizes[grepl("TN",names(sim_set_sizes))]))

#### initialize objects to start the loop and store results in later
# initialize noise to be added gradually (sd of Gaussian)
noise_seq <- seq(0, 1.5, 0.1)
# determine how many iterations will be run per scenario
n_replicates <- 10
iter <- length(noise_seq)

# create objects to store results
sim_result <- data.frame(matrix(nrow = iter, ncol = 6))
sim_result_FDR <- data.frame(matrix(nrow = iter, ncol = 6))
results_df <- data.frame()

# set colnames
colnames(sim_result) <- c("sig_hits", "overlap", "non_overlap", "sig_hits_corrected", "overlap_corrected", "non_overlap_corrected")
colnames(sim_result_FDR) <- c("sig_hits", "overlap", "non_overlap", "sig_hits_corrected", "overlap_corrected", "non_overlap_corrected")

# create objects to store enrichment results
ORA_res <- data.frame(matrix(nrow = iter, ncol=4))
colnames(ORA_res) <- c("p_val", "p_FDR", "p_val_corrected", "p_FDR_corrected")
# prepare cluster for faster computation
cl <- makeCluster(24)
clusterEvalQ(cl, library(msm))

# use a counter for quick indexing of results
counter <- 1

for (noise in noise_seq){
  set.seed(13)
  replicate_results <- data.frame()
  for (i in 1:n_replicates){
    # ASSUMPTION: CpGs within a certain set behave as they did in the set (ct is 0.1-0.9 and non-ct is < 0.1 | > 0.9)
    # select X CpGs from a set of interest (start with NeuN)
    
    # initialize matrices to hod beta values
    NeuN_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num)
    IRF8_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num)
    SOX10_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num)
    TN_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num)
    
    # simulate beta values for CpGs in the set of interest
    cpg_means <- runif(cpg_in_set_num, min = 0.15, max = 0.85)
    high_means <- runif(cpg_in_set_num, min = 0.15, max = 0.85)
    low_means <- runif(cpg_in_set_num, min = 0.01, max = 0.09)
    
    full_seq <- 1:cpg_in_set_num
    # shuffle the indices
    full_seq <- full_seq[sample(full_seq, length(full_seq))] 
    
    split_indices <- c(0, cumsum(sim_set_sizes))
    
    # Split the sequence based on the indices
    CpG_seq <- list(
      # Single sets
      NeuN = full_seq[(split_indices[1] + 1):(split_indices[2])],
      IRF8 = full_seq[(split_indices[2] + 1):(split_indices[3])],
      SOX10 = full_seq[(split_indices[3] + 1):(split_indices[4])],
      TN = full_seq[(split_indices[4] + 1):(split_indices[5])],
      
      # Sets shared between 2 cell types
      NeuN_SOX10 = full_seq[(split_indices[5] + 1):(split_indices[6])],
      NeuN_IRF8 = full_seq[(split_indices[6] + 1):(split_indices[7])],
      NeuN_TN = full_seq[(split_indices[7] + 1):(split_indices[8])],
      SOX10_IRF8 = full_seq[(split_indices[8] + 1):(split_indices[9])],
      SOX10_TN = full_seq[(split_indices[9] + 1):(split_indices[10])],
      IRF8_TN = full_seq[(split_indices[10] + 1):(split_indices[11])],
      
      # sets shared between 3 cell types
      NeuN_SOX10_IRF8 = full_seq[(split_indices[11] + 1):(split_indices[12])],
      NeuN_SOX10_TN = full_seq[(split_indices[12] + 1):(split_indices[13])],
      NeuN_IRF8_TN = full_seq[(split_indices[13] + 1):(split_indices[14])],
      SOX10_IRF8_TN = full_seq[(split_indices[14] + 1):(split_indices[15])]
    )
    
    
    # NeuN
    # get all CpG indices assigned to NeuN
    NeuN_seq <- unique(unlist(CpG_seq[grepl("NeuN", names(CpG_seq))]))
    # generate random standard deviations for each selected CpG
    sd <- runif(length(NeuN_seq), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-cpg_means[NeuN_seq])/0.5) -1)
    # fill matrix rows with truncated normal samples centered around CpG means
    NeuN_set_CpGs[NeuN_seq,] <- replicate(sample_num, rtruncnorm(length(NeuN_seq), a = 0, b = 1, mean = cpg_means[NeuN_seq], sd = sd_hetero))
    
    # get the remaining CpGs that are not part of the NeuN set
    Neun_other_seq <- full_seq[!full_seq %in% NeuN_seq]  # CpGs that are in another set
    # randomly select how many of them will be high vs. low methylated
    num <- sample(1:length(Neun_other_seq), 1)
    # randomly shuffle the remaining CpGs
    Neun_other_seq <- Neun_other_seq[sample(1:length(Neun_other_seq))]

    # split into high and low methylation groups
    NeuN_high_other_seq <- Neun_other_seq[1:num]
    NeuN_low_other_seq <- Neun_other_seq[(num+1):length(Neun_other_seq)]
    # generate high methylation means between 0.91 and 0.99
    high_means <- runif(length(NeuN_high_other_seq), min = 0.91, max = 0.99)
    # generate low methylation means between 0.01 and 0.09
    low_means <- runif(length(NeuN_low_other_seq), min = 0.01, max = 0.09)
    
    # generate variability for high means
    sd <- runif(length(high_means), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    # fill matrix for high methylation CpGs with truncated normal samples
    NeuN_set_CpGs[NeuN_high_other_seq, ] <- replicate(sample_num, rtruncnorm(length(NeuN_high_other_seq), a = 0.9, b = 1, mean = high_means, sd = sd_hetero))
    # generate variability for low means
    sd <- runif(length(low_means), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    # fill matrix for low methylation CpGs with truncated normal samples
    NeuN_set_CpGs[NeuN_low_other_seq, ] <- replicate(sample_num, rtruncnorm(length(NeuN_low_other_seq), a = 0, b = 0.1, mean = low_means, sd = sd_hetero))
    
    # -----------------------------------------------
    # Repeat the same process for IRF8, SOX10, and TN
    # -----------------------------------------------
    
    # IRF8
    IRF8_seq <- unique(unlist(CpG_seq[grepl("IRF8", names(CpG_seq))]))
    sd <- runif(length(IRF8_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[IRF8_seq])/0.5) -1)
    IRF8_set_CpGs[IRF8_seq,] <- replicate(sample_num, rtruncnorm(length(IRF8_seq), a = 0, b = 1, mean = cpg_means[IRF8_seq], sd = sd_hetero))
    
    IRF8_other_seq <- full_seq[!full_seq %in% IRF8_seq]
    num <- sample(1:length(IRF8_other_seq), 1)
    IRF8_other_seq <- IRF8_other_seq[sample(1:length(IRF8_other_seq))]
    IRF8_high_other_seq <- IRF8_other_seq[1:num]
    IRF8_low_other_seq <- IRF8_other_seq[(num+1):length(IRF8_other_seq)]
    high_means <- runif(length(IRF8_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(IRF8_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    IRF8_set_CpGs[IRF8_high_other_seq, ] <- replicate(sample_num, rtruncnorm(length(IRF8_high_other_seq), a = 0, b = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    IRF8_set_CpGs[IRF8_low_other_seq, ] <- replicate(sample_num, rtruncnorm(length(IRF8_low_other_seq), a = 0, b = 1, mean = low_means, sd = sd_hetero))
    
    # SOX10
    SOX10_seq <- unique(unlist(CpG_seq[grepl("SOX10", names(CpG_seq))]))
    sd <- runif(length(SOX10_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[SOX10_seq])/0.5) -1)
    SOX10_set_CpGs[SOX10_seq, ] <- replicate(sample_num, rtruncnorm(length(SOX10_seq), a = 0, b = 1, mean = cpg_means[SOX10_seq], sd = sd_hetero))
    
    SOX10_other_seq <- full_seq[!full_seq %in% SOX10_seq]
    num <- sample(1:length(SOX10_other_seq), 1)
    SOX10_other_seq <- SOX10_other_seq[sample(1:length(SOX10_other_seq))]
    SOX10_high_other_seq <- SOX10_other_seq[1:num]
    SOX10_low_other_seq <- SOX10_other_seq[(num+1):length(SOX10_other_seq)]
    high_means <- runif(length(SOX10_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(SOX10_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    SOX10_set_CpGs[SOX10_high_other_seq, ] <- replicate(sample_num, rtruncnorm(length(SOX10_high_other_seq), a = 0, b = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    SOX10_set_CpGs[SOX10_low_other_seq, ] <- replicate(sample_num, rtruncnorm(length(SOX10_low_other_seq), a = 0, b = 1, mean = low_means, sd = sd_hetero))
    
    # TN
    TN_seq <- unique(unlist(CpG_seq[grepl("TN", names(CpG_seq))]))
    sd <- runif(length(TN_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[TN_seq])/0.5) -1)
    TN_set_CpGs[TN_seq, ] <- replicate(sample_num, rtruncnorm(length(TN_seq), a = 0, b = 1, mean = cpg_means[TN_seq], sd = sd_hetero))
    
    TN_other_seq <- full_seq[!full_seq %in% TN_seq]
    num <- sample(1:length(TN_other_seq), 1)
    TN_other_seq <- TN_other_seq[sample(1:length(TN_other_seq))]
    TN_high_other_seq <- TN_other_seq[1:num]
    TN_low_other_seq <- TN_other_seq[(num+1):length(TN_other_seq)]
    high_means <- runif(length(TN_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(TN_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    TN_set_CpGs[TN_high_other_seq, ] <- replicate(sample_num, rtruncnorm(length(TN_high_other_seq), a = 0, b = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    TN_set_CpGs[TN_low_other_seq, ] <- replicate(sample_num, rtruncnorm(length(TN_low_other_seq), a = 0, b = 1, mean = low_means, sd = sd_hetero))
    

    #### simulate CpGs outside that set
    # create a distribution of mean beta values to sample the CpGs from
    cpg_means <- runif(cpg_out_set_num*cell_types_num, min = 0.05, max = 0.95)
    
    # generate a random sequence of means and sample from those (outside of replicate function so each CpG has the same mean across)
    rand_means <- sample(1:cpg_out_set_num)
    # generate random standard deviations for each selected CpG
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    
    clusterExport(cl, "sd_hetero")
    # create matrix of means replicated across samples for efficiency
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num), ncol = sample_num)
    # generate truncated normal beta values for each CpG and sample
    NeuN_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))

    # -----------------------------------------------
    # Repeat the same process for IRF8, SOX10, and TN
    # -----------------------------------------------
    
    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num), ncol = sample_num)
    IRF8_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))

    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num), ncol = sample_num)
    SOX10_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))

    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num), ncol = sample_num)
    TN_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))

    
    # combine set and non-set CpGs to a full betas object
    NeuN_CpGs <- rbind(NeuN_set_CpGs, NeuN_non_set_CpGs)
    IRF8_CpGs <- rbind(IRF8_set_CpGs, IRF8_non_set_CpGs)
    SOX10_CpGs <- rbind(SOX10_set_CpGs, SOX10_non_set_CpGs)
    TN_CpGs <- rbind(TN_set_CpGs, TN_non_set_CpGs)
    
    # generate cell type proportions from a Dirichlet distribution
    alpha <- rep(1, cell_types_num)
    ct_proportion <- rdirichlet(sample_num, alpha = alpha)

    
    # generate shifted cell type proportions
    alpha <- c(1,3,3,3)
    ct_proportion_change <- rdirichlet(sample_num, alpha = alpha)

    # construct the bulk values from cell type betas and proportions
    NeuN_bulk <- t(apply(NeuN_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion[,1]
    }))
    IRF8_bulk <- t(apply(IRF8_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion[,2]
    }))
    SOX10_bulk <- t(apply(SOX10_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion[,3]
    }))
    TN_bulk <- t(apply(TN_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion[,4]
    }))
    
    # combine all beta values as a linear combination of these 4 cell types
    bulk_data <- NeuN_bulk + IRF8_bulk + SOX10_bulk + TN_bulk
    rownames(bulk_data) <- paste("CpG", 1:nrow(bulk_data), sep = "_")
    
    NeuN_bulk <- t(apply(NeuN_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion_change[,1]
    }))
    IRF8_bulk <- t(apply(IRF8_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion_change[,2]
    }))
    SOX10_bulk <- t(apply(SOX10_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion_change[,3]
    }))
    TN_bulk <- t(apply(TN_CpGs, MARGIN = 1, FUN = function(x) {
      x * ct_proportion_change[,4]
    }))
    
    colMeans(ct_proportion_change)
    # save the proportions to use for correction later
    proportions_var <- rbind(ct_proportion, ct_proportion_change)
    
    # add noise to the cell type proportion estimates
    if (noise != 0){
      # generate shifted cell type proportions
      proportions_var_noise <- logit(proportions_var)
      # create a matrix with noise for each cell type value
      noise_mat <- matrix(rnorm(length(proportions_var_noise), mean = 0, sd = noise), nrow = nrow(proportions_var_noise), ncol = ncol(proportions_var_noise))
      
      # add the matrix and transform the cell type proportions back
      proportions_var_noise <- proportions_var_noise + noise_mat
      proportions_var_noise <- inverse_logit(proportions_var_noise)
      # normalize the proportions so they sum back to 1
      proportions_var_noise <- t(apply(proportions_var_noise, 1, function(x) {
        x/sum(x)
      }))
    } else{
      proportions_var_noise <- proportions_var
    }
    
    
    # combine all beta values as a linear combination of these 4 cell types
    bulk_data_change <- NeuN_bulk + IRF8_bulk + SOX10_bulk + TN_bulk
    rownames(bulk_data_change) <- paste("CpG", 1:nrow(bulk_data_change), sep = "_")
    
    # assign group labels
    case_control <- c(rep("C", sample_num), rep("AD", sample_num))
    clusterExport(cl, c("proportions_var_noise", "case_control"))
    # finally combine betas to one object to run regression on
    betas <- cbind(bulk_data, bulk_data_change)
    
    # linear regression for each CpG (row)
    lm_res <- pbapply(betas,1, FUN = lm_bulk, cl = cl, case_control)
    # correct p-values for multiple testing
    lm_res <- rbind(lm_res, p.adjust(as.numeric(lm_res[4,]), method = "BH"))
    
    set_labels <- matrix(nrow = 5, ncol = ncol(lm_res))
    set_labels[1,NeuN_seq] <- "NeuN_set"
    set_labels[2,IRF8_seq] <- "IRF8_set"
    set_labels[3,SOX10_seq] <- "SOX10_set"
    set_labels[4,TN_seq] <- "TN_set"
    set_labels[5,setdiff(1:(cpg_in_set_num + cpg_out_set_num),Reduce(union, list(NeuN_seq, IRF8_seq, SOX10_seq, TN_seq)))] <- "non_set"
    lm_res <- rbind(lm_res, set_labels)
    # assign colnames
    colnames(lm_res) <- paste("CpG", 1:ncol(lm_res), sep = "_")
    # extract significant CpGs without correction
    sig_CpGs <- (lm_res[, which(as.numeric(lm_res[4,]) < 0.05)])
    sig_CpGs_FDR <- (lm_res[, which(as.numeric(lm_res[5,]) < 0.05)])
    
    # record nr of significant hits
    sim_result[counter, "sig_hits"] <-sum(as.numeric(lm_res[4,]) < 0.05)
    sim_result_FDR[counter, "sig_hits"] <- sum(as.numeric(lm_res[5,]) < 0.05)
    
    # record nr of overlapping and non-overlapping hits
    sim_result[counter, "overlap"] <- sum(!grepl("non", sig_CpGs[10,]))
    sim_result[counter, "non_overlap"] <- sum(grepl("non", sig_CpGs[10,]))
    
    sim_result_FDR[counter, "overlap"] <- sum(!grepl("non", sig_CpGs_FDR[10,]))
    sim_result_FDR[counter, "non_overlap"] <- sum(grepl("non", sig_CpGs_FDR[10,]))
    
    ##### enrichment analysis 
    NeuN_overlap <- sum(grepl("NeuN", sig_CpGs_FDR[6,]))
    
    ORA_res[counter,"p_val"] <- phyper(q = NeuN_overlap - 1 ,
                                       m = length(sig_CpGs)/10,
                                       n = nrow(bulk_data) - length(sig_CpGs)/10,
                                       k = NeuN_set_size, lower.tail = FALSE)
    
    ORA_res[counter,"p_FDR"] <- phyper(q = NeuN_overlap - 1 ,
                                       m = length(sig_CpGs_FDR)/10,
                                       n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                       k = NeuN_set_size, lower.tail = FALSE)
    
    if (NeuN_overlap > 0) {
      
      a <- NeuN_overlap
      b <- NeuN_set_size - NeuN_overlap
      c <- ncol(sig_CpGs_FDR) - NeuN_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (NeuN_overlap == NeuN_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter, "OR"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter, "OR"] <- 0
    }
    
    IRF8_overlap <- sum(grepl("IRF8", sig_CpGs_FDR[7,]))
    
    ORA_res[counter+1,"p_val"] <- phyper(q = IRF8_overlap - 1 ,
                                         m = length(sig_CpGs)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs)/10,
                                         k = IRF8_set_size, lower.tail = FALSE)
    
    ORA_res[counter+1,"p_FDR"] <- phyper(q = IRF8_overlap - 1 ,
                                         m = length(sig_CpGs_FDR)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                         k = IRF8_set_size, lower.tail = FALSE)
    if (IRF8_overlap > 0) {
      
      a <- IRF8_overlap
      b <- IRF8_set_size - IRF8_overlap
      c <- ncol(sig_CpGs_FDR) - IRF8_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (IRF8_overlap == IRF8_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter + 1, "OR"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter + 1, "OR"] <- 0
    }
    
    SOX10_overlap <- sum(grepl("SOX10", sig_CpGs_FDR[8,]))
    
    
    ORA_res[counter+2,"p_val"] <- phyper(q = SOX10_overlap - 1 ,
                                         m = length(sig_CpGs)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs)/10,
                                         k = SOX10_set_size, lower.tail = FALSE)
    
    ORA_res[counter+2,"p_FDR"] <- phyper(q = SOX10_overlap - 1 ,
                                         m = length(sig_CpGs_FDR)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                         k = cpg_in_set_num, lower.tail = FALSE)
    if (SOX10_overlap > 0) {
      
      a <- SOX10_overlap
      b <- SOX10_set_size - SOX10_overlap
      c <- ncol(sig_CpGs_FDR) - SOX10_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (SOX10_overlap == SOX10_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter + 2, "OR"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter + 2, "OR"] <- 0
    }
    
    
    TN_overlap <- sum(grepl("TN", sig_CpGs_FDR[9,]))
    
    ORA_res[counter+3,"p_val"] <- phyper(q = TN_overlap - 1 ,
                                         m = length(sig_CpGs)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs)/10,
                                         k = TN_set_size, lower.tail = FALSE)
    
    ORA_res[counter+3,"p_FDR"] <- phyper(q = TN_overlap - 1 ,
                                         m = length(sig_CpGs_FDR)/10,
                                         n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                         k = TN_set_size, lower.tail = FALSE)
    if (TN_overlap > 0) {
      
      a <- TN_overlap
      b <- TN_set_size - TN_overlap
      c <- ncol(sig_CpGs_FDR) - TN_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (TN_overlap == TN_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter + 3, "OR"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter + 3, "OR"] <- 0
    }
    
    # correct the results for known cell type proportions
    lm_res <- pbapply(betas,1, FUN = lm_bulk_corrected, cl = cl, case_control, proportions_var_noise)
    # correct p-values for multiple testing
    lm_res <- rbind(lm_res, p.adjust(as.numeric(lm_res[4,]), method = "BH"))
    # assign set and non_set labels for ORA later
    set_labels <- matrix(nrow = 5, ncol = ncol(lm_res))
    set_labels[1,NeuN_seq] <- "NeuN_set"
    set_labels[2,IRF8_seq] <- "IRF8_set"
    set_labels[3,SOX10_seq] <- "SOX10_set"
    set_labels[4,TN_seq] <- "TN_set"
    set_labels[5,setdiff(1:(cpg_in_set_num + cpg_out_set_num),Reduce(union, list(NeuN_seq, IRF8_seq, SOX10_seq, TN_seq)))] <- "non_set"
    lm_res <- rbind(lm_res, set_labels)
    
    #lm_res <- rbind(lm_res, c(set_labels, rep("non_set", cpg_out_set_num)))
    # assign colnames
    colnames(lm_res) <- paste("CpG", 1:ncol(lm_res), sep = "_")
    # extract significant CpGs without correction
    sig_CpGs <- lm_res[, which(as.numeric(lm_res[4,]) < 0.05)]
    sig_CpGs_FDR <- lm_res[, which(as.numeric(lm_res[5,]) < 0.05)]
    
    # record nr of significant hits
    sim_result[counter, "sig_hits_corrected"] <- sum(as.numeric(lm_res[4,]) < 0.05)
    sim_result_FDR[counter, "sig_hits_corrected"] <- sum(as.numeric(lm_res[5,]) < 0.05)
    
    # record nr of overlapping and non-overlapping hits
    sim_result[counter, "overlap_corrected"] <-  sum(!grepl("non", sig_CpGs[10,]))
    sim_result[counter, "non_overlap_corrected"] <- sum(grepl("non", sig_CpGs[10,]))
    
    sim_result_FDR[counter, "overlap_corrected"] <- sum(!grepl("non", sig_CpGs_FDR[10,]))
    sim_result_FDR[counter, "non_overlap_corrected"] <- sum(grepl("non", sig_CpGs_FDR[10,]))
    
    ##### enrichment analysis 
    NeuN_overlap <- sum(grepl("NeuN", sig_CpGs_FDR[6,]))
    
    ORA_res[counter,"p_val_corrected"] <- phyper(q = NeuN_overlap - 1 ,
                                                 m = length(sig_CpGs)/10,
                                                 n = nrow(bulk_data) - length(sig_CpGs)/10,
                                                 k = NeuN_set_size, lower.tail = FALSE)
    
    ORA_res[counter,"p_FDR_corrected"] <- phyper(q = NeuN_overlap - 1 ,
                                                 m = length(sig_CpGs_FDR)/10,
                                                 n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                                 k = NeuN_set_size, lower.tail = FALSE)
    if (NeuN_overlap > 0) {
      
      a <- NeuN_overlap
      b <- NeuN_set_size - NeuN_overlap
      c <- ncol(sig_CpGs_FDR) - NeuN_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (NeuN_overlap == NeuN_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter, "OR_corrected"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter, "OR_corrected"] <- 0
    }
    
    
    IRF8_overlap <- sum(grepl("IRF8", sig_CpGs_FDR[7,]))
    
    ORA_res[counter+1,"p_val_corrected"] <- phyper(q = IRF8_overlap - 1 ,
                                                   m = length(sig_CpGs)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs)/10,
                                                   k = IRF8_set_size, lower.tail = FALSE)
    
    ORA_res[counter+1,"p_FDR_corrected"] <- phyper(q = IRF8_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                                   k = IRF8_set_size, lower.tail = FALSE)
    if (IRF8_overlap > 0) {
      # Define contingency table components
      a <- IRF8_overlap
      b <- IRF8_set_size - IRF8_overlap
      c <- (ncol(sig_CpGs_FDR)) - IRF8_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      # Apply Haldane-Anscombe correction if all test set elements are in input list
      if (IRF8_overlap == IRF8_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      # Compute corrected odds ratio
      ORA_res[counter+1, "OR_corrected"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter+1, "OR_corrected"] <- 0
    }
    
    
    SOX10_overlap <- sum(grepl("SOX10", sig_CpGs_FDR[8,]))
    
    
    ORA_res[counter+2,"p_val_corrected"] <- phyper(q = SOX10_overlap - 1 ,
                                                   m = length(sig_CpGs)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs)/10,
                                                   k = SOX10_set_size, lower.tail = FALSE)
    
    ORA_res[counter+2,"p_FDR_corrected"] <- phyper(q = SOX10_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                                   k = cpg_in_set_num, lower.tail = FALSE)
    if (SOX10_overlap > 0) {
      
      # Define contingency table components
      a <- SOX10_overlap
      b <- SOX10_set_size - SOX10_overlap
      c <- (ncol(sig_CpGs_FDR)) - SOX10_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      # Apply Haldane-Anscombe correction if all test set elements are in input list
      if (SOX10_overlap == SOX10_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      # Compute corrected odds ratio
      ORA_res[counter+2, "OR_corrected"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter+2, "OR_corrected"] <- 0
    }
    
    
    
    TN_overlap <- sum(grepl("TN", sig_CpGs_FDR[9,]))
    
    ORA_res[counter+3,"p_val_corrected"] <- phyper(q = TN_overlap - 1 ,
                                                   m = length(sig_CpGs)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs)/10,
                                                   k = TN_set_size, lower.tail = FALSE)
    
    ORA_res[counter+3,"p_FDR_corrected"] <- phyper(q = TN_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR)/10,
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR)/10,
                                                   k = TN_set_size, lower.tail = FALSE)
    if (TN_overlap > 0) {
      
      a <- TN_overlap
      b <- TN_set_size - TN_overlap
      c <- ncol(sig_CpGs_FDR) - TN_overlap
      d <- nrow(bulk_data) - (a + b + c)
      
      if (TN_overlap == TN_set_size) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }
      
      ORA_res[counter+3, "OR_corrected"] <- (a * d) / (b * c)
      
    } else {
      ORA_res[counter+3, "OR_corrected"] <- 0
    }
    
    
    # summarize all results from this iteration
    sim_result <- sim_result[complete.cases(sim_result),]
    sim_result_FDR <- sim_result_FDR[complete.cases(sim_result_FDR),]
    
    # extract results and append to the results df
    replicate_row <- cbind(sim_result, sim_result_FDR, ORA_res[1,],ORA_res[2,],ORA_res[3,],ORA_res[4,], noise)
    replicate_results <- rbind(replicate_results, replicate_row)
    
    print(paste("Computing:", noise, "at replicate number:", i))
  }
  # summarize results for all iterations
  result_summary <- data.frame(
    noise_level = noise,
    mean = apply(replicate_results, 2, mean),
    median = apply(replicate_results, 2, median),
    sd = apply(replicate_results, 2, sd),
    min = apply(replicate_results, 2, min),
    max = apply(replicate_results, 2, max),
    ci_lower = apply(replicate_results, 2, function(x) quantile(x, 0.025, na.rm = T)),
    ci_upper = apply(replicate_results, 2, function(x) quantile(x, 0.975, na.rm = T))
  )
  
  results_df <- rbind(results_df, replicate_results)
}
# stop the workers
stopCluster(cl)
# save the data
save.image("sim_workspaces/DNAm_sim_CT_lvl2.RData")