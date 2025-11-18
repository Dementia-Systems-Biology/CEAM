source("Scripts/Functions/Functions.R")
source("Scripts/Functions/plotting_functions.R")

#### Required packages
check_and_load_libraries(c(
  "ComplexUpset",
  "UpSetR",
  "patchwork",
  "purrr",
  "khroma",
  "scales",
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "furrr"
))

#### read significant CpGs from AD EWAS meta-analysis
AD_CpGs <- read.xlsx("Data/41467_2021_23243_MOESM3_ESM.xlsx", sheet = 1, startRow = 2, colNames = T)
AD_CpGs <- AD_CpGs[-1,]

#### background

# load the previously prepared 450K background
background_450 <- read.table("CpG sets/background_450k.txt", sep = "\n", header = F)
background_450 <- as.character(background_450$V1)

#### Cell type enrichment analysis

# High-specificity results
ORA_results_1 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 1)
ORA_results_1 <- ORA_results_1[!grepl("\\." ,rownames(ORA_results_1)), ] # keep only the results of the default CpG lists
upset_lvl1 <- plot_CpG_UpSet(ORA_results_1, min_set_size = 0)
# change the axis title
upset_lvl1[[3]] <- upset_lvl1[[3]] + ylab("Overlap Size")
# Medium-specificity results
ORA_results_2 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 2)
upset_lvl2 <- plot_CpG_UpSet(ORA_results_2)
# change the axis title
upset_lvl2[[3]] <- upset_lvl2[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance

# Low-specificity results
ORA_results_3 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 3)
upset_lvl3 <- plot_CpG_UpSet(ORA_results_3)
# change the axis title
upset_lvl3[[3]] <- upset_lvl3[[3]] + ylab("Overlap Size")

# use FDR to correct for multiple testing across all levels
pval_vec <- c(ORA_results_1[,"pvalue"],
              ORA_results_2[,"pvalue"],
              ORA_results_3[,"pvalue"])
FDR_vec <- p.adjust(pval_vec, method = "BH")

# append the pvalues to the ORA results
ORA_results_1$FDR <- FDR_vec[1:4]
ORA_results_2$FDR <- FDR_vec[5:8]
ORA_results_3$FDR <- FDR_vec[9:12]

# add asterisks manually to indicate significance
upset_lvl1 <- annotate_sig(upset_lvl1, sig_label = "*")
upset_lvl2 <- annotate_sig(upset_lvl2, sig_label = c("*", "", "", "*"))
upset_lvl3 <- annotate_sig(upset_lvl3, sig_label = c("*", "", "", "*"))

# combine to one plot object
combined <- ((wrap_elements(upset_lvl1)|wrap_elements(upset_lvl2)) /
 (wrap_elements(upset_lvl3) | plot_spacer()))+ 
  plot_annotation(tag_levels = 'A')

ggsave("Figures/Upsets_AD.svg", plot = combined, device = "svg", width = 160, height = 160, units = "mm")
ggsave("Figures/S6.svg", plot = combined, device = "svg", width = 160, height = 160, units = "mm")

##### permutations

# load the manifest for CpG context
manifest <- read.csv("Data/humanmethylation450_15017482_v1-2.csv", header = T, skip = 7)

levels <- c(1,2,3)
ORA_res_list <- list(ORA_results_1, ORA_results_2, ORA_results_3)

for (level in levels){
  set.seed(42)
  
  n_perm <- 1000
  background_450K <- background_450
  manifest_selected <- manifest %>%
    select(IlmnID,
           Relation_to_UCSC_CpG_Island,
           UCSC_RefGene_Group,
           Infinium_Design_Type)
  
  plan(multisession, workers = 14)
  
  perm_res <- future_map(1:n_perm, function(x) {
    rand_cpgs <- matched_sample(input_cpgs = AD_CpGs$Probe, background = background_450K, manifest = manifest_selected)
    res <- CpG_ORA(rand_cpgs, background = background_450K, specificity_level = level)
    
    data.frame(
      perm = rep(x, nrow(res)),
      OR = res$Odds.ratio
    )
  })
  
  perm_df <- bind_rows(perm_res)
  perm_df$label <- rep(rownames(ORA_res_list[[level]]),n_perm)
  
  perm_df_filt <- perm_df[!grepl("\\.", perm_df$label, ignore.case = T),]
  
  
  
  
  
  
  
  plots <- plot_OR_dist_loop(perm_df_filt, ORA_res = ORA_res_list[[level]][!grepl("\\.", rownames(ORA_res_list[[level]])),])
  
  
  names(plots)[grepl("IRF8", names(plots), ignore.case = T)] <- "MG"
  names(plots)[grepl("SOX10", names(plots), ignore.case = T)] <- "OLIG"
  names(plots)[grepl("NeuN", names(plots), ignore.case = T)] <- "NEU"
  names(plots)[grepl("TN", names(plots), ignore.case = T)] <- "AST"
  
  plots_renamed <- lapply(names(plots), function(x) {
    p <- plots[[x]]
    p + labs(title = x)
  })
  names(plots_renamed) <- names(plots)
  filename <- paste("Figures/OR_perm_AD_", level ,".RDS", sep = "")
  saveRDS(plots_renamed, file = filename)
}

##### combining all permutation plots to a supplementary figure

# loading all plots
plots1 <- readRDS("Figures/OR_perm_AD_1.rds")
plots2 <- readRDS("Figures/OR_perm_AD_2.rds")
plots3 <- readRDS("Figures/OR_perm_AD_3.rds")


combined <- 
  (plots1[["MG"]] | plots1[["OLIG"]] | plots1[["NEU"]] | plots1[["AST"]]) /
  (plots2[["MG"]] | plots2[["OLIG"]] | plots2[["NEU"]] | plots2[["AST"]]) /
  (plots3[["MG"]] | plots3[["OLIG"]] | plots3[["NEU"]] | plots3[["AST"]]) +
  plot_annotation(tag_levels = "A")


combined <- combined &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
  )

# remove the axis labels that arent necessary
for (i in seq(1:3)) {
  combined[[i]][[1]] <- combined[[i]][[1]] + theme(axis.title.y = element_text(angle = 90))
}

for (i in seq(1:4)) {
  combined[[3]][[i]] <- combined[[3]][[i]] + theme(axis.title.x = element_text())
}

ggsave("Figures/S7.svg", plot = combined, device = "svg", width = 210, height = 160, units = "mm")
