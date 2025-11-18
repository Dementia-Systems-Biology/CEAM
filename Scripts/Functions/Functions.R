# Function to check, install, and load libraries (note that this function changes the default repository temporarily)
check_and_load_libraries <- function(libraries) {
  options(repos ="https://packagemanager.posit.co/cran/2025-06-01") # an older version of the packages was required due to some unforeseen errors encountered in plotting
  
  for (lib in libraries) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      message(paste("Installing missing package:", lib))
      BiocManager::install(lib, dependencies = TRUE) # an older version of the packages was required due to some unforeseen errors encountered in plotting
    }
    library(lib, character.only = TRUE)
    message(paste("Successfully loaded package:", lib))
  }
  options(repos ="https://cran.rstudio.com/")
  
}

# function wrapper to detect hypo- and hyper-methylated CpGs with an apply function for efficiency
detect_invariant <- function(df, method, threshold){
  median_cols <- grepl("median", x = colnames(df), ignore.case = T)
  if (method == "higher"){
    return(df[apply(df[, median_cols] > threshold, MARGIN = 1, FUN = all), "cpg"])
  } else if (method == "lower"){
    return(df[apply(df[, median_cols] < threshold, MARGIN = 1, FUN = all), "cpg"])
  } else {
    stop("no method is selected")
  }
}

# function wrapper to remove hypo- and hyper-methylated CpGs with an apply function for efficiency
remove_invariant <- function(df, hyper_threshold, hypo_threshold){
  median_cols <- grepl("median", x = colnames(df), ignore.case = T)
  df <- df[!apply(df[, median_cols] > hyper_threshold, MARGIN = 1, FUN = all),]
  df <- df[!apply(df[, median_cols] < hypo_threshold, MARGIN = 1, FUN = all),]
  return(df)
}

# function to rename objects by adding a prefix
rename_obj <- function(object_name, prefix) {
  
  # Get the object
  object <- get(object_name, envir = .GlobalEnv)
  
  # Create the new name with the prefix
  new_name <- paste0(prefix, object_name)
  
  # Assign the object to the new name in the global environment
  assign(new_name, object, envir = .GlobalEnv)
  
  # Remove the original object
  rm(list = object_name, envir = .GlobalEnv)
  
  return(invisible(new_name)) # Return the new name silently
}


# function to perform linear regression on simulation data
lm_bulk <- function(row, group){
  inputdata <- data.frame(row = row, group = case_control)
  cpg_lm <- lm(formula = row ~ group, data = inputdata)
  sum_cpg <- c(summary(cpg_lm)$coefficients["groupC", ])
  return(sum_cpg)
}

# function to perform linear regression on simulation data with correction for cell type composition
lm_bulk_corrected <- function(row, group, proportions){
  inputdata <- data.frame(row = row, group = case_control,
                          NeuN = proportions[,1],
                          IRF8 = proportions[,2],
                          SOX10 = proportions[,3],
                          TN = proportions[,4])
  
  cpg_lm <- lm(formula = row ~ group + NeuN + IRF8 + SOX10 + TN, data = inputdata)
  sum_cpg <- c(summary(cpg_lm)$coefficients["groupC", ])
  return(sum_cpg)
}

logit <- function(x){
  res <- log(x/(1-x))
  return(res)
}
inverse_logit <- function(x){
  res <- exp(x)/(1+exp(x))
  return(res)
}


# function to do ORA with the custom CpG sets
CpG_ORA <- function(input, background, specificity_level) {
  
  # selects the CpG sets to be used based on the confidence level parameter
  if (specificity_level == 1){
    load("CpG sets/Lvl1_CpG_sets_Annotated.RData")
    
    signature_list <- list("IRF8_lvl1_set", "SOX10_lvl1_set", "NeuN_lvl1_set", "TN_lvl1_set")
  } else if (specificity_level == 2) {
    load("CpG sets/Lvl2_CpG_sets_Annotated.RData")
    
    signature_list <- list("IRF8_lvl2_set", "SOX10_lvl2_set", "NeuN_lvl2_set", "TN_lvl2_set")
  } else if (specificity_level == 3) {
    load("CpG sets/Lvl3_CpG_sets_Annotated.RData")
    
    signature_list <- list("IRF8_lvl3_set", "SOX10_lvl3_set", "NeuN_lvl3_set", "TN_lvl3_set")
  }
  
  
  
  
  
  if (specificity_level == 1) { # due to the subsets in high-specificity sets these need to be treated differently
    results <- data.frame(row.names = names(c(get(signature_list[[1]]),get(signature_list[[2]]),get(signature_list[[3]]),get(signature_list[[4]]))))
    for (y in signature_list){
      signature <- get(y)
      for (x in names(signature)){
        set <- get(y)[[x]]
        set <- set$cpg
        # filter the signature to be present in the background (especially important when working with non-EPIC arrays)
        set <- intersect(set, background)
        # pre-compute the intersect and union as they will be used repeatedly
        set_intersect <- intersect(input, set)
        set_union <- union(input, set)
        # find which probe identifiers overlap and drive enrichment
        results[x, "enriched_CpGs"] <- paste0(set_intersect, collapse = ",")
        # compute odds ratio
        if (((length(set) - length(set_intersect)) * (length(input) - length(set_intersect))) == 0 | length(set_intersect) == 0) { #in case the denominator is 0 or there is no overlap
          results[x, "Odds.ratio"] <- NA
        } else {
          results[x, "Odds.ratio"] <- ((length(background) - length(set_union)) * length(set_intersect))/
            ((length(set) - length(set_intersect)) * (length(input) - length(set_intersect)))
        }
        # hypergeometric test
        results[x, "pvalue"] <- phyper(q = length(set_intersect) - 1 ,
                                       m = length(input),
                                       n = length(background) - length(input), 
                                       k = length(set), lower.tail = FALSE)
        # number of CpGs overlapping
        results[x, "overlap"] <- length(intersect(input, set))
        # expected number of CpGs at random
        results[x, "expected count"] <- (length(set)/length(background))*length(input)
        # report the set size (as this may change depending on the background)
        results[x, "set_size"] <- length(set)
        
      }
      
    }
  }else{
    results <- data.frame(row.names = names(signature_list))
    for (x in signature_list){
      set <- get(x)
      set <- set$cpg
      # filter the signature to be present in the background (especially important when working with non-EPIC arrays)
      set <- intersect(set, background)
      # pre-compute the intersect and union as they will be used repeatedly
      set_intersect <- intersect(input, set)
      set_union <- union(input, set)
      
      # find which probe identifiers overlap and drive enrichment
      results[x, "enriched_CpGs"] <- paste0(set_intersect, collapse = ",")
      # compute odds ratio
      results[x, "Odds.ratio"] <- (as.numeric(length(background) - length(set_union)) * (length(set_intersect)))/
        ((length(set) - length(set_intersect)) * (length(input) - length(set_intersect)))
      
      # hypergeometric test
      results[x, "pvalue"] <- phyper(q = length(set_intersect) - 1 ,
                                     m = length(input),
                                     n = length(background) - length(input), 
                                     k = length(set), lower.tail = FALSE)
      
      # number of CpGs overlapping
      results[x, "overlap"] <- length(intersect(input, set))
      # expected number of CpGs at random
      results[x, "expected count"] <- (length(set)/length(background))*length(input)
      # report the set size (as this may change depending on the background)
      results[x, "set_size"] <- length(set)
      
    }
  }
  # adjust for multiple testing (4 tests in each specificity level)
  results$qvalue <- p.adjust(results$pvalue, "bonferroni")
  # NOTE: this assumes that a single specificity level is used to interpret the results.
  #   if multiple specificity levels will be used we implore the users to correct the raw p-values appropriately (for example BH correction)

  # returns a dataframe with the results
  return(results)
}


# function to extract relevant information from the dataframes resulting from the simulations
extract_sim_results <- function(data, group, OR_col, p_val_col){
  data %>%
    dplyr::select(all_of(c(group, OR_col, p_val_col))) %>%
        mutate(
      "{OR_col}" := ifelse(.data[[OR_col]] == 0, NA, .data[[OR_col]])
    ) %>%
  group_by(.data[[group]]) %>%
    summarise(
      median.OR = median(.data[[OR_col]], na.rm = TRUE),
      median.p = median(.data[[p_val_col]], na.rm = TRUE),
      lower_OR = quantile(.data[[OR_col]], probs = 0.25, na.rm = TRUE),
      upper_OR = quantile(.data[[OR_col]], probs = 0.75, na.rm = TRUE),
      lower_P = quantile(.data[[p_val_col]], probs = 0.25, na.rm = TRUE),
      upper_P = quantile(.data[[p_val_col]], probs = 0.75, na.rm = TRUE)
    ) %>%
    ungroup()
}


# function to resample a random background with CpGs matched by context and probe design
matched_sample <- function(input_cpgs,background, manifest) {
  # subset the manifest to only the relevant CpGs
  probe_context_bg <- manifest %>%
    filter(IlmnID %in% background)
  
  # split the multiple options for UCSC_RefGene_Group
  probe_context_bg <- probe_context_bg %>%
    tidyr::separate_rows(UCSC_RefGene_Group, sep = ";") %>%  # expands rows by each annotation
    mutate(UCSC_RefGene_Group = trimws(UCSC_RefGene_Group))    # clean up spaces
  
  # Step 3: create stratum (context grouping)
  probe_context_bg <- probe_context_bg %>%
    mutate(stratum = interaction(Relation_to_UCSC_CpG_Island, 
                                 Infinium_Design_Type, 
                                 UCSC_RefGene_Group, 
                                 drop = TRUE)) %>%
    distinct() # remove exact duplicates
  
  # get the CpG context for the real input CpGs
  context_real <- probe_context_bg %>%
    filter(IlmnID %in% input_cpgs) %>%
    slice_sample(prop = 1) %>% #randomize order so a random CpG gets removed when resolving duplicates
    distinct(IlmnID, stratum) %>%   # ensure each CpG only counted once per stratum
    count(stratum, name = "n")

  matched_samples <- context_real %>%
    group_by(stratum) %>%
    group_modify(~ {
      candidates <- probe_context_bg %>%
        filter(stratum == .y$stratum, !(IlmnID %in% input_cpgs)) %>%
        slice_sample(prop = 1) %>%       # randomize order for duplicate removal
        distinct(IlmnID, .keep_all = TRUE) # remove duplicate CpGs
      n_to_sample <- min(nrow(candidates), .x$n)
      sample_n(candidates, n_to_sample) %>%
        select(-stratum)  # <- drop grouping variable
    }) %>%
    ungroup() %>%
    pull(IlmnID)
  
  return(matched_samples)
}