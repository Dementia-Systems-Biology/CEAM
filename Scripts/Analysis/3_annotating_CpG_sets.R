source("Scripts/Functions.R")

# load required packages
check_and_load_libraries(c(
  "dplyr"
))
# load the constructed CpG sets
load("CpG_ORA_sets.RData")

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

manifest <- read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", header = T, skip = 7)

# add additional information to sets using the EPIC manifest
SOX10_lvl1_set <- lapply(X = SOX10_lvl1_set, FUN = annotate_sets, manifest = manifest)
IRF8_lvl1_set <- lapply(X = IRF8_lvl1_set, FUN = annotate_sets, manifest = manifest)
NeuN_lvl1_set <- lapply(X = NeuN_lvl1_set, FUN = annotate_sets, manifest = manifest)
TN_lvl1_set <- lapply(X = TN_lvl1_set, FUN = annotate_sets, manifest = manifest)

# save CpG sets as separate files
save(SOX10_lvl1_set, IRF8_lvl1_set, NeuN_lvl1_set, TN_lvl1_set ,file = "Lvl1_CpG_sets_Annotated.RData")

# annotate medium-specificity sets
SOX10_lvl2_set <- annotate_sets(SOX10_lvl2_set, manifest)
IRF8_lvl2_set <- annotate_sets(IRF8_lvl2_set, manifest)
NeuN_lvl2_set <- annotate_sets(NeuN_lvl2_set, manifest)
TN_lvl2_set <- annotate_sets(TN_lvl2_set, manifest)

# save CpG sets as separate files
save(SOX10_lvl2_set, IRF8_lvl2_set, NeuN_lvl2_set, TN_lvl2_set ,file = "Lvl2_CpG_sets_Annotated.RData")

# annotate low-specificity sets
SOX10_lvl3_set <- annotate_sets(SOX10_lvl3_set, manifest)
IRF8_lvl3_set <- annotate_sets(IRF8_lvl3_set, manifest)
NeuN_lvl3_set <- annotate_sets(NeuN_lvl3_set, manifest)
TN_lvl3_set <- annotate_sets(TN_lvl3_set, manifest)

# save CpG sets as separate files
save(SOX10_lvl3_set, IRF8_lvl3_set, NeuN_lvl3_set, TN_lvl3_set ,file = "Lvl3_CpG_sets_Annotated.RData")
