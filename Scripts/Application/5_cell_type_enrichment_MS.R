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

#### read data from Tiane A, Schepers M, Reijnders RA, van Veggel L, Chenine S, Rombaut B, et al. From methylation to myelination: epigenomic and transcriptomic profiling of chronic inactive demyelinated multiple sclerosis lesions. Acta Neuropathol. 2023 Aug 1;146(2):283–99. 
# significant CpGs from MS EWAS 
MS_cpgs <- read.xlsx("401_2023_2596_MOESM3_ESM.xlsx", startRow = 2, colNames = T)
# custom background from which significant CpGs were found
background <- read.csv("CpGs_NvsL.csv",header = T, row.names = 1)

#### Cell type enrichment analysis

# High-specificity results
ORA_results_1 <- CpG_ORA(input = MS_cpgs$CpG.probe, background = background$IlmnID, specificity_level = 1)

# Medium-specificity results
ORA_results_2 <- CpG_ORA(input = MS_cpgs$CpG.probe, background = background$IlmnID, specificity_level = 2)
lvl2 <- plot_CpG_UpSet(ORA_results_2, min_set_size = 150, num_breaks = 5)
# change the axis title
lvl2[[3]] <- lvl2[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
lvl2 <- annotate_sig(lvl2, sig_label = c("*", "*", "*", "*"))

# Low-specificity results
ORA_results_3 <- CpG_ORA(input = MS_cpgs$CpG.probe, background = background$IlmnID, specificity_level = 3)
lvl3 <- plot_CpG_UpSet(ORA_results_3, min_set_size = 150, num_breaks = 5)
# change the axis title
lvl3[[3]] <- lvl3[[3]] + ylab("Overlap Size")
# add asterisks manually to indicate significance
lvl3 <- annotate_sig(lvl3, sig_label = c("*", "*", "*", "*"))

# combine to one plot object
combined <- (wrap_elements(lvl2) | wrap_elements(lvl3)) + 
  plot_annotation(tag_levels = "A")