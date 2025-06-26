source("Scripts/Functions.R")
# load required packages
check_and_load_libraries(c(
  "dplyr",
  "pbapply",
  "parallel"
))
# load the betas and summary statistics for both cohorts (ensure they are differently named)
load("UKBBN2_betas.RData")
load("UKBBN2_summary_results.RData")
summary_res_UKBBN2 <- summary_res
rm(summary_res)

load("BDR_betas.RData")
load("BDR_summary_results.RData")
summary_res_BDR <- summary_res
rm(summary_res)

# this function finds intermediately methylated median betas in each cell type
find_intermediate_CpGs <- function(df, celltype_col, other_col, threshold_range, min_diff, multiple_cell_cpg){
  #' @param df A data.frame containing summary statistics of methylation data per column
  #' @param celltype_col the column name containing the summary statistic of the target cell type
  #' @param other_col a vector of column names containing the summary statistics of the other cell types
  #' @param threshold_range a single number specifying the width of the central threshold (eg. 0.4 means between 0.1 and 0.9)
  #' @param min_diff minimum difference between target cell type methylation and all others for inclusion
  #' @param multiple_cell_cpg logical indicating if cell type-specific CpGs should be included for multiple cell types
  

  # CpGs that are interemdiately methylated and less methylated than the other 3 cell types
  idx <- parApply(cl = cl, df[, other_col] > 0.5 + threshold_range & 
                 abs(0.5 - df[, celltype_col]) < threshold_range & 
                 df[, celltype_col] < df[,other_col] - min_diff, 1, all)
  
  # CpGs that are interemdiately methylated and more methylated than the other 3 cell types
  idx2 <- parApply(cl = cl, df[, other_col] < 0.5 - threshold_range & 
                  abs(0.5 - df[, celltype_col]) < threshold_range & 
                  df[, celltype_col] > df[,other_col] + min_diff, 1, all)
  
  # CpGs that are interemdiately methylated and in between the methylation values of other cell types
  idx9 <- parApply(cl = cl, abs(0.5 - df[,other_col]) > threshold_range & 
                  abs(0.5 - df[, celltype_col]) < threshold_range & 
                  abs(df[, celltype_col] - df[,other_col]) > min_diff, 1, all)
  
  if (multiple_cell_cpg == T){
    # finding intermediately methylated CpGs in the other cell types
    idx3 <- abs(0.5 - df[,other_col]) < threshold_range
    idx3 <- cbind(idx3, abs(0.5 - df[,celltype_col]) < threshold_range)
    
    # find which rows have the target cell type and all others except one within the threshold range
    selected_rows <- rowSums(idx3[, -ncol(idx3)] > 0) & idx3[, ncol(idx3)] # splits it into all cols except the last one & only the last column
    
    # set all rows where these conditions are not fulfilled to entirely FALSE
    idx3[!selected_rows, ] <- FALSE 
    
    # find which cell type is outside the threshold range
    idx4 <- !idx3
    
    # set all rows where these conditions aren not fulfilled to entirely FALSE
    idx4[!selected_rows,] <- FALSE 
    
    # subset the dataframe to compare values more easily
    df_subset <- df[,c(other_col, celltype_col)]
    
    # set all values that are FALSE to 0
    values_idx3 <- df_subset * idx3
    values_idx4 <- df_subset * idx4
    values_idx4[values_idx4 == 0] <- NA
    
    # find the extreme values idx4 either contains 1 or 2 values since 3 or 2 cells are non-extreme
    extreme_values <- t(parApply(cl = cl,values_idx4, 1, range, na.rm = T))
    extreme_values <- ifelse(is.finite(extreme_values), extreme_values, NA) # set all the infinite to NA so they dont get computed with the diff later
   
     # compute the difference between each non-zero value of the celltypes within threshold range
    diff_matrix1 <- abs(values_idx3 - extreme_values[,1])
    diff_matrix2 <- abs(values_idx3 - extreme_values[,2])
    
    # fix subtracting the 0 value in idx3 from its complementary extreme value in idx4
    diff_matrix1[diff_matrix1 == extreme_values[,1]] <- NA # removes the extreme value being subtracted by 0 in each row
    diff_matrix2[diff_matrix2 == extreme_values[,2]] <- NA # removes the extreme value being subtracted by 0 in each row
    
    # find where all celltypes in threshold range have at least the min_diff with the hyper- or hypomethylated celltype
    idx5 <- parApply(cl = cl, diff_matrix1, 1, function(x) {
      if (all(is.na(x))) {
        return(FALSE)  # If all values in the row are NA, return FALSE
      }
      all(x > min_diff, na.rm = TRUE)  # Otherwise, check if any value is > 0.1
    })
    idx6 <- parApply(cl = cl, diff_matrix2, 1, function(x) {
      if (all(is.na(x))) {
        return(FALSE)  # If all values in the row are NA, return FALSE
      }
      all(x > min_diff, na.rm = TRUE)  # Otherwise, check if any value is > 0.1
    })
    idx7 <- idx5 & idx6
  } else{
    idx7 <- NULL
  }
  
  
  return(list(low_cpgs = df[idx,], high_cpgs = df[idx2, ], multi_cell_cpgs = df[idx7, ], flip_flop = df[idx9,]))
}

##### High-Specificity #####
# parallelize the operations for slight increase in speed
cl <- makeCluster(14)
clusterExport(cl, "summary_res_UKBBN2")
SOX10_UKBBN2_lvl1 <- find_intermediate_CpGs(summary_res_UKBBN2, celltype_col = "Sox10_Median", other_col = c("TN_Median","IRF8_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
IRF8_UKBBN2_lvl1 <- find_intermediate_CpGs(summary_res_UKBBN2, celltype_col = "IRF8_Median", other_col = c("TN_Median","Sox10_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
NeuN_UKBBN2_lvl1 <- find_intermediate_CpGs(summary_res_UKBBN2, celltype_col = "NeuN_Median", other_col = c("TN_Median","IRF8_Median","Sox10_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
TN_UKBBN2_lvl1 <- find_intermediate_CpGs(summary_res_UKBBN2, celltype_col = "TN_Median", other_col = c("Sox10_Median","IRF8_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)

clusterExport(cl, "summary_res_BDR")
SOX10_BDR_lvl1 <- find_intermediate_CpGs(summary_res_BDR, celltype_col = "Sox10_Median", other_col = c("TN_Median","IRF8_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
IRF8_BDR_lvl1 <- find_intermediate_CpGs(summary_res_BDR, celltype_col = "IRF8_Median", other_col = c("TN_Median","Sox10_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
NeuN_BDR_lvl1 <- find_intermediate_CpGs(summary_res_BDR, celltype_col = "NeuN_Median", other_col = c("TN_Median","IRF8_Median","Sox10_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)
TN_BDR_lvl1 <- find_intermediate_CpGs(summary_res_BDR, celltype_col = "TN_Median", other_col = c("Sox10_Median","IRF8_Median","NeuN_Median"), threshold_range = 0.4, min_diff = 0.1, multiple_cell_cpg = T)


# extract the CpG names of all the CpGs that are found in both cohorts with the same patterns (in both cohorts low, high or in between)
SOX10_lvl1_set <- list(SOX10_lvl1.1 = intersect(SOX10_UKBBN2_lvl1$low_cpgs$cpg, SOX10_BDR_lvl1$low_cpgs$cpg),
                       SOX10_lvl1.9 = intersect(SOX10_UKBBN2_lvl1$high_cpgs$cpg, SOX10_BDR_lvl1$high_cpgs$cpg),
                       SOX10_lvl1.5 = setdiff(intersect(SOX10_UKBBN2_lvl1$flip_flop$cpg, SOX10_BDR_lvl1$flip_flop$cpg),
                                              c(intersect(SOX10_UKBBN2_lvl1$low_cpgs$cpg, SOX10_BDR_lvl1$low_cpgs$cpg),
                                                intersect(SOX10_UKBBN2_lvl1$high_cpgs$cpg, SOX10_BDR_lvl1$high_cpgs$cpg))),
                       SOX10_lvl1 = unique(c(intersect(SOX10_UKBBN2_lvl1$low_cpgs$cpg, SOX10_BDR_lvl1$low_cpgs$cpg),
                                             intersect(SOX10_UKBBN2_lvl1$high_cpgs$cpg, SOX10_BDR_lvl1$high_cpgs$cpg),
                                             intersect(SOX10_UKBBN2_lvl1$flip_flop$cpg, SOX10_BDR_lvl1$flip_flop$cpg))))

IRF8_lvl1_set <- list(IRF8_lvl1.1 = intersect(IRF8_UKBBN2_lvl1$low_cpgs$cpg, IRF8_BDR_lvl1$low_cpgs$cpg),
                      IRF8_lvl1.9 = intersect(IRF8_UKBBN2_lvl1$high_cpgs$cpg, IRF8_BDR_lvl1$high_cpgs$cpg),
                      IRF8_lvl1.5 = setdiff(intersect(IRF8_UKBBN2_lvl1$flip_flop$cpg, IRF8_BDR_lvl1$flip_flop$cpg),
                                            c(intersect(IRF8_UKBBN2_lvl1$low_cpgs$cpg, IRF8_BDR_lvl1$low_cpgs$cpg),
                                              intersect(IRF8_UKBBN2_lvl1$high_cpgs$cpg, IRF8_BDR_lvl1$high_cpgs$cpg))),
                      IRF8_lvl1 = unique(c(intersect(IRF8_UKBBN2_lvl1$low_cpgs$cpg, IRF8_BDR_lvl1$low_cpgs$cpg),
                                           intersect(IRF8_UKBBN2_lvl1$high_cpgs$cpg, IRF8_BDR_lvl1$high_cpgs$cpg),
                                           intersect(IRF8_UKBBN2_lvl1$flip_flop$cpg, IRF8_BDR_lvl1$flip_flop$cpg))))

NeuN_lvl1_set <- list(NeuN_lvl1.1 = intersect(NeuN_UKBBN2_lvl1$low_cpgs$cpg, NeuN_BDR_lvl1$low_cpgs$cpg),
                      NeuN_lvl1.9 = intersect(NeuN_UKBBN2_lvl1$high_cpgs$cpg, NeuN_BDR_lvl1$high_cpgs$cpg),
                      NeuN_lvl1.5 = setdiff(intersect(NeuN_UKBBN2_lvl1$flip_flop$cpg, NeuN_BDR_lvl1$flip_flop$cpg),
                                            c(intersect(NeuN_UKBBN2_lvl1$low_cpgs$cpg, NeuN_BDR_lvl1$low_cpgs$cpg),
                                              intersect(NeuN_UKBBN2_lvl1$high_cpgs$cpg, NeuN_BDR_lvl1$high_cpgs$cpg))),
                      NeuN_lvl1 = unique(c(intersect(NeuN_UKBBN2_lvl1$low_cpgs$cpg, NeuN_BDR_lvl1$low_cpgs$cpg),
                                           intersect(NeuN_UKBBN2_lvl1$high_cpgs$cpg, NeuN_BDR_lvl1$high_cpgs$cpg),
                                           intersect(NeuN_UKBBN2_lvl1$flip_flop$cpg, NeuN_BDR_lvl1$flip_flop$cpg))))

TN_lvl1_set <- list(TN_lvl1.1 = intersect(TN_UKBBN2_lvl1$low_cpgs$cpg, TN_BDR_lvl1$low_cpgs$cpg),
                    TN_lvl1.9 = intersect(TN_UKBBN2_lvl1$high_cpgs$cpg, TN_BDR_lvl1$high_cpgs$cpg),
                    TN_lvl1.5 = setdiff(intersect(TN_UKBBN2_lvl1$flip_flop$cpg, TN_BDR_lvl1$flip_flop$cpg),
                                        c(intersect(TN_UKBBN2_lvl1$low_cpgs$cpg, TN_BDR_lvl1$low_cpgs$cpg),
                                          intersect(TN_UKBBN2_lvl1$high_cpgs$cpg, TN_BDR_lvl1$high_cpgs$cpg))),
                    TN_lvl1 = unique(c(intersect(TN_UKBBN2_lvl1$low_cpgs$cpg, TN_BDR_lvl1$low_cpgs$cpg),
                                       intersect(TN_UKBBN2_lvl1$high_cpgs$cpg, TN_BDR_lvl1$high_cpgs$cpg),
                                       intersect(TN_UKBBN2_lvl1$flip_flop$cpg, TN_BDR_lvl1$flip_flop$cpg))))

##### Medium-Specificity #####

# check the shared cpg overlap (has to overlap with the same 2 and 3 cell types in both cohorts)
shared_cpg_UKBBN2 <- data.frame(cpg = summary_res_UKBBN2$cpg)
shared_cpg_UKBBN2$SOX10 <- ifelse(shared_cpg_UKBBN2$cpg %in% SOX10_UKBBN2_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_UKBBN2$IRF8 <- ifelse(shared_cpg_UKBBN2$cpg %in% IRF8_UKBBN2_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_UKBBN2$NeuN <- ifelse(shared_cpg_UKBBN2$cpg %in% NeuN_UKBBN2_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_UKBBN2$TN <- ifelse(shared_cpg_UKBBN2$cpg %in% TN_UKBBN2_lvl1$multi_cell_cpgs$cpg , 1, 0)

shared_cpg_BDR <- data.frame(cpg = summary_res_BDR$cpg)
shared_cpg_BDR$SOX10 <- ifelse(shared_cpg_BDR$cpg %in% SOX10_BDR_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_BDR$IRF8 <- ifelse(shared_cpg_BDR$cpg %in% IRF8_BDR_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_BDR$NeuN <- ifelse(shared_cpg_BDR$cpg %in% NeuN_BDR_lvl1$multi_cell_cpgs$cpg , 1, 0)
shared_cpg_BDR$TN <- ifelse(shared_cpg_BDR$cpg %in% TN_BDR_lvl1$multi_cell_cpgs$cpg , 1, 0)

shared_cpg_joined <- inner_join(shared_cpg_UKBBN2, shared_cpg_BDR, by = "cpg")
shared_cpg <- filter(shared_cpg_joined, SOX10.x == SOX10.y, IRF8.x == IRF8.y, NeuN.x == NeuN.y, TN.x == TN.y)
shared_cpg <- shared_cpg[rowSums(shared_cpg[,-1]) > 0, ]
shared_cpg_names <- shared_cpg$cpg

# extract the same results but adding the CpGs shared between 3 cell types
SOX10_lvl2_set <- unique(c(SOX10_lvl1_set$SOX10_lvl1,
                           intersect(intersect(SOX10_UKBBN2_lvl1$multi_cell_cpgs$cpg, SOX10_BDR_lvl1$multi_cell_cpgs$cpg), shared_cpg_names)))

IRF8_lvl2_set <- unique(c(IRF8_lvl1_set$IRF8_lvl1,
                          intersect(intersect(IRF8_UKBBN2_lvl1$multi_cell_cpgs$cpg, IRF8_BDR_lvl1$multi_cell_cpgs$cpg), shared_cpg_names)))

NeuN_lvl2_set <- unique(c(NeuN_lvl1_set$NeuN_lvl1,
                          intersect(intersect(NeuN_UKBBN2_lvl1$multi_cell_cpgs$cpg, NeuN_BDR_lvl1$multi_cell_cpgs$cpg), shared_cpg_names)))

TN_lvl2_set <- unique(c(TN_lvl1_set$TN_lvl1,
                        intersect(intersect(TN_UKBBN2_lvl1$multi_cell_cpgs$cpg, TN_BDR_lvl1$multi_cell_cpgs$cpg), shared_cpg_names)))


##### Low-Specificity #####
generate_lvl_3_sets <- function(df, median_col, FDR_cols, other_median_col, threshold_range, p_val, min_diff){
  #' @param df A data.frame containing summary statistics of methylation data per column
  #' @param median_col the column name containing the summary statistic of the target cell type
  #' @param FDR_cols the columns containing the FDR-corrected p-values for each pairwise comparison
  #' @param other_col a vector of column names containing the summary statistics of the other cell types
  #' @param threshold_range a single number specifying the width of the central threshold (eg. 0.4 means between 0.1 and 0.9)
  #' @param p_val a single number specifying p-value threshold considered significant
  #' @param min_diff minimum difference between target cell type methylation and all others for inclusion

  idx <- parApply(cl = cl, df[, FDR_cols] < p_val & # significant difference between target cell type and all others
                 abs(0.5 - df[, median_col]) < threshold_range & # intermediately methylated
                 abs(df[, median_col] - df[,other_median_col]) > min_diff, 1, all) # at least 0.1 difference between target cell type and all others
  return(df[idx,])
}

SOX10_UKBBN2_lvl3 <- generate_lvl_3_sets(summary_res_UKBBN2, 
                                         median_col = "Sox10_Median", 
                                         FDR_cols = c("Sox10_TN_FDR_p-value", "Sox10_IRF8_FDR_p-value", "Sox10_NeuN_FDR_p-value"), 
                                         other_median_col = c("TN_Median", "IRF8_Median", "NeuN_Median"), 
                                         threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

IRF8_UKBBN2_lvl3 <- generate_lvl_3_sets(summary_res_UKBBN2, 
                                        median_col = "IRF8_Median", 
                                        FDR_cols = c("TN_IRF8_FDR_p-value", "Sox10_IRF8_FDR_p-value", "IRF8_NeuN_FDR_p-value"), 
                                        other_median_col = c("TN_Median", "Sox10_Median", "NeuN_Median"), 
                                        threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

NeuN_UKBBN2_lvl3 <- generate_lvl_3_sets(summary_res_UKBBN2, 
                                        median_col = "NeuN_Median", 
                                        FDR_cols = c("IRF8_NeuN_FDR_p-value", "TN_NeuN_FDR_p-value", "Sox10_NeuN_FDR_p-value"), 
                                        other_median_col = c("TN_Median", "IRF8_Median", "Sox10_Median"), 
                                        threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

TN_UKBBN2_lvl3 <- generate_lvl_3_sets(summary_res_UKBBN2, 
                                      median_col = "TN_Median", 
                                      FDR_cols = c("Sox10_TN_FDR_p-value", "TN_IRF8_FDR_p-value", "TN_NeuN_FDR_p-value"), 
                                      other_median_col = c("Sox10_Median", "IRF8_Median", "NeuN_Median"), 
                                      threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)
# extract for BDR
SOX10_BDR_lvl3 <- generate_lvl_3_sets(summary_res_BDR, 
                                      median_col = "Sox10_Median", 
                                      FDR_cols = c("Sox10_IRF8_FDR_p-value", "Sox10_NeuN_FDR_p-value", "Sox10_TN_FDR_p-value"), 
                                      other_median_col = c("TN_Median", "IRF8_Median", "NeuN_Median"), 
                                      threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

IRF8_BDR_lvl3 <- generate_lvl_3_sets(summary_res_BDR, 
                                     median_col = "IRF8_Median", 
                                     FDR_cols = c("Sox10_IRF8_FDR_p-value", "IRF8_NeuN_FDR_p-value", "IRF8_TN_FDR_p-value"), 
                                     other_median_col = c("TN_Median", "Sox10_Median", "NeuN_Median"), 
                                     threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

NeuN_BDR_lvl3 <- generate_lvl_3_sets(summary_res_BDR, 
                                     median_col = "NeuN_Median", 
                                     FDR_cols = c("IRF8_NeuN_FDR_p-value", "NeuN_TN_FDR_p-value", "Sox10_NeuN_FDR_p-value"), 
                                     other_median_col = c("TN_Median", "IRF8_Median", "Sox10_Median"), 
                                     threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

TN_BDR_lvl3 <- generate_lvl_3_sets(summary_res_BDR, 
                                   median_col = "TN_Median", 
                                   FDR_cols = c("Sox10_TN_FDR_p-value", "IRF8_TN_FDR_p-value", "NeuN_TN_FDR_p-value"), 
                                   other_median_col = c("Sox10_Median", "IRF8_Median", "NeuN_Median"), 
                                   threshold_range = 0.4, p_val = 0.05, min_diff = 0.1)

# keep only the CpGs found in both cohorts
SOX10_lvl3_set <- intersect(SOX10_UKBBN2_lvl3$cpg, SOX10_BDR_lvl3$cpg)
IRF8_lvl3_set <- intersect(IRF8_UKBBN2_lvl3$cpg, IRF8_BDR_lvl3$cpg)
NeuN_lvl3_set <- intersect(NeuN_UKBBN2_lvl3$cpg, NeuN_BDR_lvl3$cpg)
TN_lvl3_set <- intersect(TN_UKBBN2_lvl3$cpg, TN_BDR_lvl3$cpg)

#### combining all sets to be cumulative (so lvl3 contains at least all in lvl2 and lvl1 etc.)

# lvl2 sets should fully include lvl1 sets
length(intersect(SOX10_lvl1_set$SOX10_lvl1, SOX10_lvl2_set))/length(SOX10_lvl1_set$SOX10_lvl1) == 1
length(intersect(IRF8_lvl1_set$IRF8_lvl1, IRF8_lvl2_set)) / length(IRF8_lvl1_set$IRF8_lvl1) == 1
length(intersect(NeuN_lvl1_set$NeuN_lvl1, NeuN_lvl2_set)) / length(NeuN_lvl1_set$NeuN_lvl1) == 1
length(intersect(TN_lvl1_set$TN_lvl1, TN_lvl2_set)) / length(TN_lvl1_set$TN_lvl1) == 1

# merge lvl2 and 3
SOX10_lvl3_set <- unique(c(SOX10_lvl2_set, SOX10_lvl3_set))
IRF8_lvl3_set <- unique(c(IRF8_lvl2_set, IRF8_lvl3_set))
NeuN_lvl3_set <- unique(c(NeuN_lvl2_set, NeuN_lvl3_set))
TN_lvl3_set <- unique(c(TN_lvl2_set, TN_lvl3_set))

# lvl3 sets should then fullyinclude lvl2 sets
length(intersect(SOX10_lvl3_set, SOX10_lvl2_set))/length(SOX10_lvl2_set) == 1
length(intersect(IRF8_lvl3_set, IRF8_lvl2_set)) / length(IRF8_lvl2_set) == 1
length(intersect(NeuN_lvl3_set, NeuN_lvl2_set)) / length(NeuN_lvl2_set) == 1
length(intersect(TN_lvl3_set, TN_lvl2_set)) / length(TN_lvl2_set) == 1

save(SOX10_lvl1_set, IRF8_lvl1_set, NeuN_lvl1_set, TN_lvl1_set, 
     SOX10_lvl2_set, IRF8_lvl2_set, NeuN_lvl2_set, TN_lvl2_set,
     SOX10_lvl3_set, IRF8_lvl3_set, NeuN_lvl3_set, TN_lvl3_set, file = "CpG_ORA_sets.RData")
