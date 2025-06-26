source("Scripts/Functions.R")

check_and_load_libraries(c(
  "dplyr",
  "ggplot2",
  "pbapply",
  "parallel",
  "nlme",
  "lmerTest",
  "performance",
  "data.table"
))

load("UKBBN2_betas.RData")

# prepare the clusters
cl <- makeCluster(14) # 16 cores is fast enough, but more cores would improve speed further if desired

UKBBN2_pheno$Cell_Type <- gsub(" ", replacement = "", x = UKBBN2_pheno$Cell_Type)

cell_types <- unique(UKBBN2_pheno$Cell_Type)
pairwise_combinations <- combn(cell_types, 2)
pairwise_names <- apply(pairwise_combinations,2,paste,collapse='_')

results <- rep(list(NULL), 6)

names(results) <- pairwise_names

#### Linear Mixed Effect (LME) Models
# iterate over each pairwise comparison
for (i in 1:ncol(pairwise_combinations)){
  
  # select the two cell types to be compared
  ct1 <- pairwise_combinations[1,i]
  ct2 <- pairwise_combinations[2,i]
  
  # subset the data to obtain pheno and betas for only these cell type samples
  pheno_subset <- UKBBN2_pheno[UKBBN2_pheno$Cell_Type == ct1 | UKBBN2_pheno$Cell_Type == ct2, ]
  beta_subset <- UKBBN2_betas[ ,UKBBN2_pheno$Cell_Type == ct1 | UKBBN2_pheno$Cell_Type == ct2]
  
  # ensure that meta data and betas are in the same order
  if (all(rownames(pheno_subset) == colnames(beta_subset))){
    # create a variable that compares one cell type to all other three
    pheno_subset$ct1_vs_ct2 = ifelse(pheno_subset$Cell_Type == ct1, 1, 0)
    
    
    # extract all relevant variables for the lme function
    Cell_type = pheno_subset$ct1_vs_ct2
    individual = pheno_subset$Individual_ID
    age = pheno_subset$Age
    sex = pheno_subset$Sex
    
    # function to compute LME models and extract their results
    f_lme <- function(row, Cell_type, individual, sex, age){ # for each CpG it tests for cell type, and corrects for individual, sex and age
      
      inputdata  = data.frame(row = row, Cell_type = Cell_type, individual = individual, age = age, sex = sex)
      cpg_lme = nlme::lme(row~ age + sex + Cell_type, random=~1|individual,  control = nlme::lmeControl(opt = "optim"))
      
      # extract the relevant results
      sum_cpg = c(summary(cpg_lme)$tTable["Cell_type", ])
      return(sum_cpg)
      
    }
    # provide all necessary data (inlcuding the function) to the workers
    clusterExport(cl, c("f_lme", "Cell_type", "individual", "age", "sex"))
    
    
    # run the function for each row of the betas matrix
    AllLME = pbapply(beta_subset,1,f_lme, Cell_type, individual, age, sex, cl = cl)
    
    results[[pairwise_names[i]]] <- AllLME
  } else{
    stop("betas and pheno are not in the same order!")
  }
}

save.image("Workspaces/pairwise_LME_UKBBN2_results.RData")


summary_res <- data.frame(matrix(NA, nrow = ncol(results[[1]]), ncol = (6*4) + (4*3))) # ncol is calculated based on 3 summary statistics per cell type (4*3) and 4 summary statistics for each pairwise LME (6*4)
rownames(summary_res) <- colnames(results[[1]])

counter <- 1

for (ct in cell_types){
  beta_subset <- UKBBN2_betas[ ,UKBBN2_pheno$Cell_Type == ct]
  if (all(rownames(beta_subset) == rownames(summary_res))){
    summary_res[, counter] <- pbapply(beta_subset, cl = cl, MARGIN = 1, FUN = median)
    summary_res[, counter+1] <- rowMeans(beta_subset)
    summary_res[, counter +2] <- pbapply(beta_subset, cl = cl, MARGIN = 1, FUN = sd)
  }
  colnames(summary_res)[counter:(counter+2)] <- paste(ct, c("Median","Mean", "sd"), sep = "_") 
  
  counter <- counter + 3
}
# add a row index and CpG identifier column for easier access later
summary_res$Row.Index <- seq_len(nrow(summary_res))
summary_res$cpg <- rownames(summary_res)

# (optional) investigate how many and which CpGs are hypo- or hyper-methylated in all cell types
hypermethylated_cpg <- detect_invariant(summary_res, method = "higher", threshold = 0.95)
hypomethylated_cpg <- detect_invariant(summary_res, method = "lower", threshold = 0.05)

# remove CpGs hypo- or hyper-methylated in all cell types
summary_res <- remove_invariant(summary_res, hyper_threshold = 0.95, hypo_threshold = 0.05) # removes 229,178 CpGs

# extract statistics from LMEs
for (i in 1:length(results)) {
  
  # append the LME statistics to the matching CpGs in the summary_res
  summary_res[,counter+0] <- results[[i]]["Value", summary_res$cpg]
  summary_res[,counter+1] <- results[[i]]["Std.Error", summary_res$cpg]
  summary_res[,counter+2] <- results[[i]]["p-value", summary_res$cpg]
  summary_res[,counter+3] <- p.adjust(summary_res[,counter+2], method = "BH")
  
  # name the columns by the pairwise comparison and summary statistic
  colnames(summary_res)[counter:(counter+3)] <- paste(names(results)[i], c("Value", "std.Error", "p-value", "FDR_p-value"), 
                                                      sep = "_")
  # increase the counter to add the next 4 columns
  counter <- counter + 4
  
}
# stop the workers
stopCluster(cl)

save(summary_res, file = "UKBBN2_summary_results.RData")
