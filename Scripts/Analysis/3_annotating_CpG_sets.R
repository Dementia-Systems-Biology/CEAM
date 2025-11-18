source("Scripts/Functions/Functions.R")

# load required packages
check_and_load_libraries(c(
  "dplyr",
  "qs"
))
# load the constructed CpG sets
load("Data/CpG_ORA_sets.RData")

annotate_sets <- function(set, manifest){
  results <- data.frame(cpg = set)
  # select relevant information from manifest
  selected_cols <- manifest[, c("Name", "UCSC_RefGene_Name", 
                                "UCSC_RefGene_Group", 
                                "UCSC_CpG_Islands_Name",
                                "Relation_to_UCSC_CpG_Island", 
                                "Regulatory_Feature_Group")]
  # join the CpG set and manifest by their probe identifiers
  left_join(results, y = selected_cols, by = join_by(cpg == Name))
}

manifest <- read.csv("Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv", header = T, skip = 7)

# add additional information to sets using the EPIC manifest
SOX10_lvl1_set <- lapply(X = SOX10_lvl1_set, FUN = annotate_sets, manifest = manifest)
IRF8_lvl1_set <- lapply(X = IRF8_lvl1_set, FUN = annotate_sets, manifest = manifest)
NeuN_lvl1_set <- lapply(X = NeuN_lvl1_set, FUN = annotate_sets, manifest = manifest)
TN_lvl1_set <- lapply(X = TN_lvl1_set, FUN = annotate_sets, manifest = manifest)

# save CpG sets as separate files
save(SOX10_lvl1_set, IRF8_lvl1_set, NeuN_lvl1_set, TN_lvl1_set ,file = "CpG sets/Lvl1_CpG_sets_Annotated.RData")

write.table(SOX10_lvl1_set$SOX10_lvl1$cpg, file = "CpG sets/OLIG_high.txt", col.names = F, row.names = F)
write.table(IRF8_lvl1_set$IRF8_lvl1$cpg, file = "CpG sets/MG_high.txt", col.names = F, row.names = F)
write.table(NeuN_lvl1_set$NeuN_lvl1$cpg, file = "CpG sets/NEU_high.txt", col.names = F, row.names = F)
write.table(TN_lvl1_set$TN_lvl1$cpg, file = "CpG sets/TN_high.txt", col.names = F, row.names = F)

# annotate medium-specificity sets
SOX10_lvl2_set <- annotate_sets(SOX10_lvl2_set, manifest)
IRF8_lvl2_set <- annotate_sets(IRF8_lvl2_set, manifest)
NeuN_lvl2_set <- annotate_sets(NeuN_lvl2_set, manifest)
TN_lvl2_set <- annotate_sets(TN_lvl2_set, manifest)

# save CpG sets as separate files
save(SOX10_lvl2_set, IRF8_lvl2_set, NeuN_lvl2_set, TN_lvl2_set ,file = "CpG sets/Lvl2_CpG_sets_Annotated.RData")
write.table(SOX10_lvl2_set$cpg, file = "CpG sets/OLIG_med.txt", col.names = F, row.names = F)
write.table(IRF8_lvl2_set$cpg, file = "CpG sets/MG_med.txt", col.names = F, row.names = F)
write.table(NeuN_lvl2_set$cpg, file = "CpG sets/NEU_med.txt", col.names = F, row.names = F)
write.table(TN_lvl2_set$cpg, file = "CpG sets/TN_med.txt", col.names = F, row.names = F)

# annotate low-specificity sets
SOX10_lvl3_set <- annotate_sets(SOX10_lvl3_set, manifest)
IRF8_lvl3_set <- annotate_sets(IRF8_lvl3_set, manifest)
NeuN_lvl3_set <- annotate_sets(NeuN_lvl3_set, manifest)
TN_lvl3_set <- annotate_sets(TN_lvl3_set, manifest)

# save CpG sets as separate files
save(SOX10_lvl3_set, IRF8_lvl3_set, NeuN_lvl3_set, TN_lvl3_set ,file = "CpG sets/Lvl3_CpG_sets_Annotated.RData")
write.table(SOX10_lvl3_set$cpg, file = "CpG sets/OLIG_low.txt", col.names = F, row.names = F)
write.table(IRF8_lvl3_set$cpg, file = "CpG sets/MG_low.txt", col.names = F, row.names = F)
write.table(NeuN_lvl3_set$cpg, file = "CpG sets/NEU_low.txt", col.names = F, row.names = F)
write.table(TN_lvl3_set$cpg, file = "CpG sets/TN_low.txt", col.names = F, row.names = F)

##### write the CpG sets to a format usable by the shiny app
cpg_list <- list(lvl1 = NA, lvl1.5 = NA, lvl2 = NA, lvl3 = NA)

# high-specificity sets (lvl1 sets)
cpg_list$lvl1 <- list(
  MG = IRF8_lvl1_set$IRF8_lvl1$cpg,
  NEU = NeuN_lvl1_set$NeuN_lvl1$cpg,
  OLIG = SOX10_lvl1_set$SOX10_lvl1$cpg,
  AST = TN_lvl1_set$TN_lvl1$cpg
)

# high-specificity sets with sub-levels (lvl1.5 sets)
cpg_list$lvl1.5 <- list(
  MG = IRF8_lvl1_set$IRF8_lvl1$cpg,
  MG_hyper = IRF8_lvl1_set$IRF8_lvl1.9$cpg,
  MG_hypo = IRF8_lvl1_set$IRF8_lvl1.1$cpg,
  
  NEU = NeuN_lvl1_set$NeuN_lvl1$cpg,
  NEU_hyper = NeuN_lvl1_set$NeuN_lvl1.9$cpg,
  NEU_hypo  = NeuN_lvl1_set$NeuN_lvl1.1$cpg,
  
  OLIG = SOX10_lvl1_set$SOX10_lvl1$cpg,
  OLIG_hyper = SOX10_lvl1_set$SOX10_lvl1.9$cpg,
  OLIG_hypo  = SOX10_lvl1_set$SOX10_lvl1.1$cpg,
  
  AST = TN_lvl1_set$TN_lvl1$cpg,
  AST_hyper = TN_lvl1_set$TN_lvl1.9$cpg,
  AST_hypo  = TN_lvl1_set$TN_lvl1.1$cpg
  )

# medium-specificity sets (lvl2 sets)
cpg_list$lvl2 <- list(
  MG = IRF8_lvl2_set$cpg,
  NEU = NeuN_lvl2_set$cpg,
  OLIG = SOX10_lvl2_set$cpg,
  AST = TN_lvl2_set$cpg
)

# low-specificity sets (lvl3 sets)
cpg_list$lvl3 <- list(
  MG = IRF8_lvl3_set$cpg,
  NEU = NeuN_lvl3_set$cpg,
  OLIG = SOX10_lvl3_set$cpg,
  AST = TN_lvl3_set$cpg
)

# save the lists of CpGs to use for the app
qsave(cpg_list, "CpG sets/CpG_sets.qs") # to inspect the CpG lists in detail we refer to the .txt files or the RData files containing the dataframes with annotation
