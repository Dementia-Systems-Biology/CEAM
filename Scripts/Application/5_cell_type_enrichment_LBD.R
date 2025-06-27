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
  "ggplot2"
))

#### read summary statistics from LB disease EWAS
load("LBMeta_PrimarySumStats.Rdata")

# select significant CpGs with suggestive P-value
sig_CpGs <- rownames(resMetaCC[resMetaCC$P_Fixed < 1e-5,])

#### Cell type enrichment analysis

# High-specificity results
ORA_results_1 <- CpG_ORA(input = sig_CpGs, background = rownames(resMetaCC), specificity_level = 1)

# Medium-specificity results
ORA_results_2 <- CpG_ORA(input = sig_CpGs, background = rownames(resMetaCC), specificity_level = 2)
lvl2 <- plot_CpG_UpSet(ORA_results_2)
# change the axis title
lvl2[[3]] <- lvl2[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
lvl2 <- annotate_sig(lvl2, sig_label = c("*", "", "", ""))


# Low-specificity results
ORA_results_3 <- CpG_ORA(input = sig_CpGs, background = rownames(resMetaCC), specificity_level = 3)
lvl3 <- plot_CpG_UpSet(ORA_results_3)
# change the axis title
lvl3[[3]] <- lvl3[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
lvl3 <- annotate_sig(lvl3, sig_label = c("*", "", "", ""))

# combine to one plot object
combined <- (wrap_elements(lvl2) | wrap_elements(lvl3)) + 
  plot_annotation(tag_levels = "A")

