source("Scripts/Functions.R")
#### Required packages
check_and_load_libraries(c(
  "ComplexUpset",
  "UpSetR",
  "patchwork",
  "purrr",
  "khroma",
  "scales",
  "dplyr",
  "ggplot2",
  "openxlsx"
))

#### read significant CpGs from AD EWAS meta-analysis
AD_CpGs <- read.xlsx("41467_2021_23243_MOESM3_ESM.xlsx", sheet = 1, startRow = 2, colNames = T)
AD_CpGs <- AD_CpGs[-1,]

# load 450K background (London1 cohort of Smith et al.)
load("background_450K.RData")

#### Cell type enrichment analysis

# High-specificity results
ORA_results_1 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 1)

# Medium-specificity results
ORA_results_2 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 2)
upset_lvl2 <- plot_CpG_UpSet(ORA_results_2, "Enrichment results lvl2")
# change the axis title
upset_lvl2[[3]] <- upset_lvl2[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
upset_lvl2 <- annotate_sig(upset_lvl2, sig_label = c("*", "", "", ""))

# Low-specificity results
ORA_results_3 <- CpG_ORA(input = AD_CpGs$Probe, background = background_450, specificity_level = 3)
upset_lvl3 <- plot_CpG_UpSet(ORA_results_3, "Enrichment results lvl3")
# change the axis title
upset_lvl3[[3]] <- upset_lvl3[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
upset_lvl3 <- annotate_sig(upset_lvl3, sig_label = c("*", "", "", "*"))

# combine to one plot object
combined <- (wrap_elements(upset_lvl2)|wrap_elements(upset_lvl3)) + 
  plot_annotation(tag_levels = 'A')