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

# find set sizes of real data
set_sizes <- c(23746, 19143, 2636, 3001) # set sizes of high-specificity sets
names(set_sizes) <- c("NeuN", "IRF8", "SOX10", "TN")

background_size <- 785614

# initialize the nr of CpGs and samples that should be simulated
sample_num <- 250
cell_types_num <- 4
scaling_factor <- 0.25

# determine the simualted set sizes based on real set sizes
NeuN_set_size <- round(scaling_factor * set_sizes[["NeuN"]])
IRF8_set_size <-  round(scaling_factor * set_sizes[["IRF8"]])
SOX10_set_size <-  round(scaling_factor * set_sizes[["SOX10"]])
TN_set_size <-  round(scaling_factor * set_sizes[["TN"]])
all_cpg_num <- background_size * scaling_factor

# find how many CpGs are in a set and how many are not
cpg_in_set_num <- NeuN_set_size + IRF8_set_size + SOX10_set_size + TN_set_size
cpg_out_set_num <- all_cpg_num - (cpg_in_set_num + IRF8_set_size + SOX10_set_size + TN_set_size)

#### initialize objects to start the loop and store results
# initialize DMPs to be added gradually (sd of Gaussian)
DMP_seq <- c(0, seq(0, 1.5, 0.1))

# determine how many iterations will be run per scenario
n_replicates <- 10
iter <- length(DMP_seq)

# create objects to store results
sim_result <- data.frame(matrix(nrow = iter, ncol = 6))
sim_result_FDR <- data.frame(matrix(nrow = iter, ncol = 6))
results_df <- data.frame()

# set colnames
colnames(sim_result) <- c("sig_hits_corrected", "overlap_corrected", "non_overlap_corrected")
colnames(sim_result_FDR) <- c("sig_hits_corrected", "overlap_corrected", "non_overlap_corrected")

# create objects to store enrichment results
ORA_res <- data.frame(matrix(nrow = iter, ncol=4))
colnames(ORA_res) <- c("p_val_corrected", "p_FDR_corrected")

#prepare cluster for faster computation
cl <- makeCluster(24)
clusterEvalQ(cl, library(msm))  # Load msm on each worker

# use a counter for quick indexing of results
counter <- 1
for (DMPs in DMP_seq){
  set.seed(13)
  
  # initialize number of DMPS for each cell type and background
  n_dmps <- DMPs
  dmp_ct <- list( # current implementation can not handle 0 in any cell type
    NeuN = 0.5,
    IRF8 = 0.20,
    SOX10 = 0.001,
    TN = 0.001,
    background = 0.1
  )
  
  
  initial_dmps <- round(n_dmps * unlist(dmp_ct))
  # Cap each DMP count at max available CpGs in that group
  dmp_alloc <- pmin(initial_dmps, list(NeuN_set_size, IRF8_set_size, SOX10_set_size, TN_set_size, cpg_out_set_num))
  remaining_dmps <- n_dmps - sum(unlist(dmp_alloc))
  # assign remaining DMPs to background
  bg_remaining_capacity <- cpg_out_set_num - dmp_alloc[["background"]]
  bg_extra_dmps <- min(remaining_dmps, bg_remaining_capacity)
  
  dmp_alloc[["background"]] <- dmp_alloc[["background"]] + bg_extra_dmps
  # initialize object to store results
  replicate_results <- data.frame()
  for (i in 1:n_replicates){
    # ASSUMPTION: CpGs within a certain set behave as they did in the set (ct is 0.1-0.9 and non-ct is < 0.1 | > 0.9)
    # select X CpGs from a set of interest (start with NeuN)
    
    # initialize matrices to hold beta values
    NeuN_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num*2)
    IRF8_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num*2)
    SOX10_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num*2)
    TN_set_CpGs <- matrix(nrow = cpg_in_set_num, ncol = sample_num*2)
    
    # simulate beta values for CpGs in the set of interest
    cpg_means <- runif(cpg_in_set_num, min = 0.15, max = 0.85)
    # change means for DMPs
    cpg_means_dmp <- inverse_logit(logit(cpg_means)  + sample(c(seq(-1,-0.5,0.05), seq(0.5,1,0.05)), length(cpg_means), replace = T))
    
    # create indices for all CpGs
    full_seq <- 1:cpg_in_set_num
    full_seq <- full_seq[sample(full_seq, length(full_seq))]
    
    # NeuN
    # select a subset of CpG sites for the NeuN-specific set
    NeuN_seq <- full_seq[1:NeuN_set_size]
    # generate random standard deviations for each selected CpG
    sd <- runif(length(NeuN_seq), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-cpg_means[NeuN_seq])/0.5) -1)
    # fill matrix rows with truncated normal samples centered around CpG means
    NeuN_set_CpGs[NeuN_seq,1:sample_num] <- replicate(sample_num, rtnorm(NeuN_set_size, lower = 0, upper = 1, mean = cpg_means[NeuN_seq], sd = sd_hetero))
    # duplicate first set of samples
    NeuN_set_CpGs[NeuN_seq,(sample_num+1):(sample_num*2)] <-  NeuN_set_CpGs[NeuN_seq,1:sample_num]
    
    # randomly select CpGs from the set to be differentially methylated 
    NeuN_dmp_seq <- sample(NeuN_seq, dmp_alloc[["NeuN"]], replace = F)
    # generate random standard deviations for each selected CpG
    sd <- runif(length(NeuN_dmp_seq), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-cpg_means[NeuN_dmp_seq])/0.5) -1)
    # simulate new values for DMPs in the second set of samples
    NeuN_set_CpGs[NeuN_dmp_seq,(sample_num+1):(sample_num*2)] <- replicate(sample_num, rtnorm(dmp_alloc[["NeuN"]], lower = 0, upper = 1, mean = cpg_means_dmp[NeuN_dmp_seq], sd = sd_hetero))
    
    
    # get the remaining CpGs that are not part of the NeuN set
    Neun_other_seq <- full_seq[!full_seq %in% NeuN_seq]  # CpGs that are in another set
    # randomly select how many of them will be high vs. low methylated
    num <- sample(1:length(Neun_other_seq)-1, 1)
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
    # simulate high methylation values
    NeuN_set_CpGs[NeuN_high_other_seq, ] <- replicate(sample_num*2, rtnorm(length(NeuN_high_other_seq), lower = 0.9, upper = 1, mean = high_means, sd = sd_hetero))
   
    # generate variability for high means
    sd <- runif(length(low_means), min = 0, max = 0.1)
    # make variance heteroskedastic
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    # simulate low methylation values
    NeuN_set_CpGs[NeuN_low_other_seq, ] <- replicate(sample_num*2, rtnorm(length(NeuN_low_other_seq), lower = 0, upper = 0.1, mean = low_means, sd = sd_hetero))
    
    # -----------------------------------------------
    # Repeat the same process for IRF8, SOX10, and TN
    # -----------------------------------------------

    # IRF8
    # all CpGs belonging to IRF8
    IRF8_seq <- full_seq[(NeuN_set_size + 1):(NeuN_set_size + IRF8_set_size)]
    sd <- runif(length(IRF8_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[IRF8_seq])/0.5) -1)
    IRF8_set_CpGs[IRF8_seq,] <- replicate(sample_num, rtnorm(IRF8_set_size, lower = 0, upper = 1, mean = cpg_means[IRF8_seq], sd = sd_hetero))
    
    # DMP CpGs belonging to IRF8
    
    # only run DMP part when DMPs are assigned to IRF8
    if (dmp_alloc[["IRF8"]] > 0){
      IRF8_dmp_seq <- sample(IRF8_seq, dmp_alloc[["IRF8"]], replace = F)
      sd <- runif(length(IRF8_dmp_seq), min = 0, max = 0.1)
      sd_hetero <- sd*abs((abs(0.5-cpg_means[IRF8_dmp_seq])/0.5) -1)
      IRF8_set_CpGs[IRF8_dmp_seq,(sample_num+1):(sample_num*2)] <- replicate(sample_num, rtnorm(dmp_alloc[["IRF8"]], lower = 0, upper = 1, mean = cpg_means_dmp[IRF8_dmp_seq], sd = sd_hetero))
    }
    # CpGs not in IRF8
    IRF8_other_seq <- full_seq[!full_seq %in% IRF8_seq]
    num <- sample(1:length(IRF8_other_seq)-1, 1)
    IRF8_other_seq <- IRF8_other_seq[sample(1:length(IRF8_other_seq))]
    IRF8_high_other_seq <- IRF8_other_seq[1:num]
    IRF8_low_other_seq <- IRF8_other_seq[(num+1):length(IRF8_other_seq)]
    high_means <- runif(length(IRF8_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(IRF8_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    IRF8_set_CpGs[IRF8_high_other_seq, ] <- replicate(sample_num*2, rtnorm(length(IRF8_high_other_seq), lower = 0, upper = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    IRF8_set_CpGs[IRF8_low_other_seq, ] <- replicate(sample_num*2, rtnorm(length(IRF8_low_other_seq), lower = 0, upper = 1, mean = low_means, sd = sd_hetero))
    
    # SOX10
    SOX10_seq <- full_seq[(NeuN_set_size + IRF8_set_size + 1):(NeuN_set_size + IRF8_set_size + SOX10_set_size)]
    sd <- runif(length(SOX10_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[SOX10_seq])/0.5) -1)
    SOX10_set_CpGs[SOX10_seq, ] <- replicate(sample_num, rtnorm(SOX10_set_size, lower = 0, upper = 1, mean = cpg_means[SOX10_seq], sd = sd_hetero))
    
    # only run DMP part when DMPs are assigned to SOX10
    if (dmp_alloc[["SOX10"]] > 0){
      SOX10_dmp_seq <- sample(SOX10_seq, dmp_alloc[["SOX10"]], replace = F)
      sd <- runif(length(SOX10_dmp_seq), min = 0, max = 0.1)
      sd_hetero <- sd*abs((abs(0.5-cpg_means[SOX10_dmp_seq])/0.5) -1)
      SOX10_set_CpGs[SOX10_dmp_seq,(sample_num+1):(sample_num*2)] <- replicate(sample_num, rtnorm(dmp_alloc[["SOX10"]], lower = 0, upper = 1, mean = cpg_means_dmp[SOX10_dmp_seq], sd = sd_hetero))
      
    }
    
    SOX10_other_seq <- full_seq[!full_seq %in% SOX10_seq]
    num <- sample(1:length(SOX10_other_seq)-1, 1)
    SOX10_other_seq <- SOX10_other_seq[sample(1:length(SOX10_other_seq))]
    SOX10_high_other_seq <- SOX10_other_seq[1:num]
    SOX10_low_other_seq <- SOX10_other_seq[(num+1):length(SOX10_other_seq)]
    high_means <- runif(length(SOX10_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(SOX10_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    SOX10_set_CpGs[SOX10_high_other_seq, ] <- replicate(sample_num*2, rtnorm(length(SOX10_high_other_seq), lower = 0, upper = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    SOX10_set_CpGs[SOX10_low_other_seq, ] <- replicate(sample_num*2, rtnorm(length(SOX10_low_other_seq), lower = 0, upper = 1, mean = low_means, sd = sd_hetero))
    
    # TN
    TN_seq <- full_seq[(NeuN_set_size + IRF8_set_size + SOX10_set_size + 1):cpg_in_set_num]
    sd <- runif(length(TN_seq), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[TN_seq])/0.5) -1)
    TN_set_CpGs[TN_seq, ] <- replicate(sample_num, rtnorm(TN_set_size, lower = 0, upper = 1, mean = cpg_means[TN_seq], sd = sd_hetero))
    
    # only run DMP part when DMPs are assigned to TN
    if (dmp_alloc[["TN"]] > 0){
      TN_dmp_seq <- sample(TN_seq, dmp_alloc[["TN"]], replace = F)
      sd <- runif(length(TN_dmp_seq), min = 0, max = 0.1)
      sd_hetero <- sd*abs((abs(0.5-cpg_means[TN_dmp_seq])/0.5) -1)
      TN_set_CpGs[TN_dmp_seq,(sample_num+1):(sample_num*2)] <- replicate(sample_num, rtnorm(dmp_alloc[["TN"]], lower = 0, upper = 1, mean = cpg_means_dmp[TN_dmp_seq], sd = sd_hetero))
    }
    
    TN_other_seq <- full_seq[!full_seq %in% TN_seq]
    num <- sample(1:length(TN_other_seq)-1, 1)
    TN_other_seq <- TN_other_seq[sample(1:length(TN_other_seq))]
    TN_high_other_seq <- TN_other_seq[1:num]
    TN_low_other_seq <- TN_other_seq[(num+1):length(TN_other_seq)]
    high_means <- runif(length(TN_high_other_seq), min = 0.91, max = 0.99)
    low_means <- runif(length(TN_low_other_seq), min = 0.01, max = 0.09)
    sd <- runif(length(high_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-high_means)/0.5) -1)
    TN_set_CpGs[TN_high_other_seq, ] <- replicate(sample_num*2, rtnorm(length(TN_high_other_seq), lower = 0, upper = 1, mean = high_means, sd = sd_hetero))
    sd <- runif(length(low_means), min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-low_means)/0.5) -1)
    TN_set_CpGs[TN_low_other_seq, ] <- replicate(sample_num*2, rtnorm(length(TN_low_other_seq), lower = 0, upper = 1, mean = low_means, sd = sd_hetero))
    

    
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
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num*2), ncol = sample_num*2)
    # generate truncated normal beta values for each CpG and sample
    NeuN_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))
    
    # -----------------------------------------------
    # Repeat the same process for IRF8, SOX10, and TN
    # -----------------------------------------------
    
    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num*2), ncol = sample_num*2)
    IRF8_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))
    

    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num*2), ncol = sample_num*2)
    SOX10_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))
    
    
    rand_means <- sample(1:cpg_out_set_num)
    sd <- runif(cpg_out_set_num, min = 0, max = 0.1)
    sd_hetero <- sd*abs((abs(0.5-cpg_means[rand_means])/0.5) -1)
    clusterExport(cl, "sd_hetero")
    mean_mat <- matrix(rep(cpg_means[rand_means], sample_num*2), ncol = sample_num*2)
    TN_non_set_CpGs <- parApply(mean_mat, 2, cl = cl, function(m) rtnorm(length(m), mean = m, sd = sd_hetero, lower = 0, upper = 1))
    
    
    # combine set and non-set CpGs to a full betas object
    NeuN_CpGs <- rbind(NeuN_set_CpGs, NeuN_non_set_CpGs)
    IRF8_CpGs <- rbind(IRF8_set_CpGs, IRF8_non_set_CpGs)
    SOX10_CpGs <- rbind(SOX10_set_CpGs, SOX10_non_set_CpGs)
    TN_CpGs <- rbind(TN_set_CpGs, TN_non_set_CpGs)
    

    # generate cell type proportions from a Dirichlet distribution
    alpha <- rep(1, cell_types_num)
    ct_proportion <- rdirichlet(sample_num*2, alpha = alpha)
    
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
    
    # add the background DMPs
    
    bg_dmp_seq <- sample((cpg_in_set_num+1):all_cpg_num, dmp_alloc[["background"]], replace = F)
    bulk_data[bg_dmp_seq, (sample_num+1):(sample_num*2)] <- inverse_logit(logit(
      bulk_data[bg_dmp_seq, 1:sample_num])  + sample(c(-1,1), length(bg_dmp_seq), replace = T))
    
    # assign group labels
    case_control <- c(rep("C", sample_num), rep("AD", sample_num))
    clusterExport(cl, c("ct_proportion", "case_control"))

    # correct the results for known cell type proportions
    lm_res <- pbapply(bulk_data,1, FUN = lm_bulk_corrected, cl = cl, case_control, ct_proportion)
    # correct p-values for multiple testing
    lm_res <- rbind(lm_res, p.adjust(as.numeric(lm_res[4,]), method = "BH"))
    # assign set and non_set labels for ORA later
    set_labels <- matrix(nrow = 1, ncol = cpg_in_set_num)
    set_labels[,NeuN_seq] <- "NeuN_set"
    set_labels[,IRF8_seq] <- "IRF8_set"
    set_labels[,SOX10_seq] <- "SOX10_set"
    set_labels[,TN_seq] <- "TN_set"
    
    lm_res <- rbind(lm_res, c(set_labels, rep("non_set", cpg_out_set_num)))
    # assign colnames
    colnames(lm_res) <- paste("CpG", 1:ncol(lm_res), sep = "_")
    # extract significant CpGs without correction
    sig_CpGs <- as.data.frame(lm_res[, which(as.numeric(lm_res[4,]) < 0.05)])
    sig_CpGs_FDR <- as.data.frame(lm_res[, which(as.numeric(lm_res[5,]) < 0.05)])
    
    # record nr of significant hits
    sim_result[counter, "sig_hits_corrected"] <- sum(as.numeric(lm_res[4,]) < 0.05)
    sim_result_FDR[counter, "sig_hits_corrected"] <- sum(as.numeric(lm_res[5,]) < 0.05)
    
    # record nr of overlapping and non-overlapping hits
    sim_result[counter, "overlap_corrected"] <-  sum(!grepl("non", sig_CpGs[6,]))
    sim_result[counter, "non_overlap_corrected"] <- sum(grepl("non", sig_CpGs[6,]))
    
    sim_result_FDR[counter, "overlap_corrected"] <- sum(!grepl("non", sig_CpGs_FDR[6,]))
    sim_result_FDR[counter, "non_overlap_corrected"] <- sum(grepl("non", sig_CpGs_FDR[6,]))
    
    ##### enrichment analysis 
    NeuN_overlap <- sum(grepl("NeuN", sig_CpGs_FDR[6,]))
    
    ORA_res[counter,"p_val_corrected"] <- phyper(q = NeuN_overlap - 1 ,
                                                 m = length(sig_CpGs),
                                                 n = nrow(bulk_data) - length(sig_CpGs),
                                                 k = NeuN_set_size, lower.tail = FALSE)
    
    ORA_res[counter,"p_FDR_corrected"] <- phyper(q = NeuN_overlap - 1 ,
                                                 m = length(sig_CpGs_FDR),
                                                 n = nrow(bulk_data) - length(sig_CpGs_FDR),
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
    
    
    IRF8_overlap <- sum(grepl("IRF8", sig_CpGs_FDR[6,]))
    
    ORA_res[counter+1,"p_val_corrected"] <- phyper(q = IRF8_overlap - 1 ,
                                                   m = length(sig_CpGs),
                                                   n = nrow(bulk_data) - length(sig_CpGs),
                                                   k = IRF8_set_size, lower.tail = FALSE)
    
    ORA_res[counter+1,"p_FDR_corrected"] <- phyper(q = IRF8_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR),
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR),
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
    
    
    SOX10_overlap <- sum(grepl("SOX10", sig_CpGs_FDR[6,]))
    
    
    ORA_res[counter+2,"p_val_corrected"] <- phyper(q = SOX10_overlap - 1 ,
                                                   m = length(sig_CpGs),
                                                   n = nrow(bulk_data) - length(sig_CpGs),
                                                   k = SOX10_set_size, lower.tail = FALSE)
    
    ORA_res[counter+2,"p_FDR_corrected"] <- phyper(q = SOX10_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR),
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR),
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
    
    
    
    TN_overlap <- sum(grepl("TN", sig_CpGs_FDR[6,]))
    
    ORA_res[counter+3,"p_val_corrected"] <- phyper(q = TN_overlap - 1 ,
                                                   m = length(sig_CpGs),
                                                   n = nrow(bulk_data) - length(sig_CpGs),
                                                   k = TN_set_size, lower.tail = FALSE)
    
    ORA_res[counter+3,"p_FDR_corrected"] <- phyper(q = TN_overlap - 1 ,
                                                   m = length(sig_CpGs_FDR),
                                                   n = nrow(bulk_data) - length(sig_CpGs_FDR),
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
    # extract results and append to the results df
    replicate_row <- cbind(sim_result, sim_result_FDR, ORA_res[1,],ORA_res[2,],ORA_res[3,],ORA_res[4,], DMPs)
    replicate_results <- rbind(replicate_results, replicate_row)
    print(paste("Computing:", DMPs, "at replicate number:", i))
    
  }
  # summarize results for all iterations
  result_summary <- data.frame(
    DMPs_level = DMPs,
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
save.image("sim_workspaces/DNAm_sim_DMP_lvl1.RData")