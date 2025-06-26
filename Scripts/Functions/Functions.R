# Function to check, install, and load libraries
check_and_load_libraries <- function(libraries) {
  for (lib in libraries) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      message(paste("Installing missing package:", lib))
      BiocManager::install(lib, dependencies = TRUE)
    }
    library(lib, character.only = TRUE)
    message(paste("Successfully loaded package:", lib))
  }
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
    load("D:/University Data/Internship_raw_data_analysis/Lvl1_CpG_sets_Annotated.RData")
    
    signature_list <- list("IRF8_lvl1_set", "SOX10_lvl1_set", "NeuN_lvl1_set", "TN_lvl1_set")
  } else if (specificity_level == 2) {
    load("D:/University Data/Internship_raw_data_analysis/Lvl2_CpG_sets_Annotated.RData")
    
    signature_list <- list("IRF8_lvl2_set", "SOX10_lvl2_set", "NeuN_lvl2_set", "TN_lvl2_set")
  } else if (specificity_level == 3) {
    load("D:/University Data/Internship_raw_data_analysis/Lvl3_CpG_sets_Annotated.RData")
    
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
        results[x, "Odds.ratio"] <- ((length(background) - length(set_union)) * length(set_intersect))/
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
  # returns a dataframe with the results
  return(results)
}

# function to create UpSet plots from cell type enrichment results
plot_CpG_UpSet <- function(ORA_result, min_set_size = 3, num_breaks = 3){
  # Prepare input list
  tmp <- strsplit(ORA_result$enriched_CpGs, split = ",")
  # change the names to the abbreviations
  names(tmp) <- sapply(rownames(ORA_result), function(x) {
    if (grepl("NeuN", x)) {
      "NEU"
    } else if (grepl("IRF8", x)) {
      "MG"
    } else if (grepl("SOX10", x)) {
      "OLIG"
    } else if (grepl("TN", x)) {
      "AST"
    } else {
      x  # fallback to original name if no match
    }
  })
  
  list_input <- UpSetR::fromList(tmp)
  colnames(list_input) <- names(tmp)
  rownames(list_input) <- unique(unlist(tmp))
  list_input <- as.data.frame(list_input == 1)
  
  # Plot
  ComplexUpset::upset(
    list_input,
    names(tmp),
    name = "CpG sets",
    base_annotations=list(
      'Intersection size'=intersection_size(counts=FALSE) # change to TRUE if values should be displayed
    ),
    min_size = min_set_size,
    set_sizes = (
      upset_set_size() +
        scale_y_reverse(
          breaks = scales::pretty_breaks(n = num_breaks) # makes nicer breaks for visualizing total overlap size
        ) +
        theme(axis.text.x = element_text(angle = 90))
    ),
    width_ratio = 0.3)
}

# function to manually annotate UpSet plots
annotate_sig <- function(upset, sig_label) {
  #' @param upset the UpSet plot object from ComplexUpset to be annotated
  #' @param sig_label a vector of characters used to indicate significance in the barplot, commonly "*"
  pbuilt <- ggplot_build(upset[[3]])
  bar_data <- pbuilt$data[[2]]  # Usually the first layer is the bars
  upset[[3]] <- upset[[3]] +  annotate("text", x = as.numeric(bar_data[,"x"]) -0.15, y = bar_data[,"count"] + 0.1*bar_data[,"count"], label = sig_label, size = 6)
  return(upset)
}