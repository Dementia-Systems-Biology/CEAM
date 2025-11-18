source("Scripts/Functions/Functions.R")
source("Scripts/Functions/plotting_functions.R")
check_and_load_libraries(c(
  "sva",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "cowplot"
))

load("Data/UKBBN2_combat_vs_nocombat.RData")
load("Data/BDR_betas.RData")

# change the way that cell type is noted in UKBBN2 pheno
UKBBN2_pheno <- UKBBN2_pheno %>%
  mutate(
    `cell type:ch1` = recode(
      `cell type:ch1`,
      "Sox10 " = "Oligodendrocytes",
      "TN" = "Astrocyte-enriched",
      "NeuN" = "Neurons",
      "IRF8" = "Microglia"
    )
  )
BDR_pheno <- BDR_pheno %>%
  mutate(
    `cell type:ch1` = recode(
      `cell type:ch1`,
      "Oligodendrocytes" = "Oligodendrocytes",
      "Astrocytes" = "Astrocyte-enriched",
      "Neurons" = "Neurons",
      "Microglia" = "Microglia"
    )
  )
# select colors (colorblind friendly)
IBM_cols <- c("#D81B60","#004D40", "#FF9307","#1E88E5")

#### PCA of BDR

# Performing PCA
pca_result <-  prcomp(t(BDR_betas),        
                      retx = TRUE,
                      center =TRUE,
                      scale = FALSE,
                      rank. = 5)
scores_only_all <- as.data.frame(pca_result$x)

#PCA results are combined with phenotype data
pca_result_sample_info <- scores_only_all %>%
  mutate(ID = row.names(scores_only_all)) %>%
  inner_join(BDR_pheno, by = c("ID" = "Basename"))

(PCA_celltype_BDR <- ggplot(pca_result_sample_info, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = `cell type:ch1` ), alpha = 0.7) +
    coord_fixed(1) +
    labs(
      x = paste0("PC1 (",round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 2), "%)"),
      color = "Cell Type") + 
    scale_color_manual(
      values = c(
        "Oligodendrocytes" = IBM_cols[1],
        "Astrocyte-enriched" = IBM_cols[4],
        "Neurons" = IBM_cols[2],
        "Microglia" = IBM_cols[3]
      )
    ))

#### PCA of UKBBN2 before ComBat
# first remove all the problematic CpGs from this as well

### Removing SNP CpG sites
SNPall = read.table("Data/SNPProbes_McCartney.txt", header = T) # obtained from McCartney DL, Walker RM, Morris SW, McIntosh AM, Porteous DJ, Evans KL. Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip. Genom Data. 2016 Sep;9:22–4. 
SNPallEURAF = SNPall[which(SNPall$EUR_AF >= 0.05 & SNPall$EUR_AF <= 0.95),]

#### Removing Cross-hybridising sites
crosshyb = read.table("Data/CrossHydridisingProbes_McCartney.txt") # obtained from McCartney DL, Walker RM, Morris SW, McIntosh AM, Porteous DJ, Evans KL. Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip. Genom Data. 2016 Sep;9:22–4. 

#### Removing X and Y chromosomes probes
manifest_EPIC = read.csv("Data/MethylationEPIC_v-1-0_B4.csv", header = T, skip = 7)
manifest_EPIC_XY = manifest_EPIC[manifest_EPIC$CHR %in% c("X", "Y"),]

to_remove_cpgs = unique(c(SNPallEURAF$IlmnID, manifest_EPIC_XY$IlmnID, crosshyb$V1))
#### From the 3 steps above you should remove 72240 CpGs in total

#### Removing Non-CpG sites
non_CpG_probes = rownames(UKBBN2_betas_nocombat[grep("ch.", rownames(UKBBN2_betas_nocombat)),])
UKBBN2_betas_nocombat = UKBBN2_betas_nocombat[!rownames(UKBBN2_betas_nocombat) %in% non_CpG_probes,]
dim(UKBBN2_betas_nocombat)

#### Removing Reference SNPs
ref_SNP_probes = rownames(UKBBN2_betas_nocombat[grep("rs", rownames(UKBBN2_betas_nocombat)),])
UKBBN2_betas_nocombat = UKBBN2_betas_nocombat[!rownames(UKBBN2_betas_nocombat) %in% to_remove_cpgs,]
UKBBN2_betas_nocombat = UKBBN2_betas_nocombat[!rownames(UKBBN2_betas_nocombat) %in% ref_SNP_probes,]

#### Removing Non-CpG sites
non_CpG_probes = rownames(UKBBN2_betas_nocombat[grep("ch.", rownames(UKBBN2_betas_nocombat)),])
UKBBN2_betas_nocombat = UKBBN2_betas_nocombat[!rownames(UKBBN2_betas_nocombat) %in% non_CpG_probes,]
dim(UKBBN2_betas_nocombat)

# Performing PCA
pca_result <-  prcomp(t(UKBBN2_betas_nocombat),        
                      retx = TRUE,
                      center =TRUE,
                      scale = FALSE,
                      rank. = 5)
scores_only_all <- as.data.frame(pca_result$x)

#PCA results are combined with phenotype data
pca_result_sample_info <- scores_only_all %>%
  mutate(ID = row.names(scores_only_all)) %>%
  inner_join(UKBBN2_pheno, by = c("ID" = "Basename"))

(PCA_celltype_UKBBN2_NC <- ggplot(pca_result_sample_info, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = `cell type:ch1` ), alpha = 0.7) +
    coord_fixed(1) +
    labs(
      x = paste0("PC1 (",round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 2), "%)"),
      color = "Cell Type") + 
    scale_color_manual(
      values = c(
        "Oligodendrocytes" = IBM_cols[1],
        "Astrocyte-enriched" = IBM_cols[4],
        "Neurons" = IBM_cols[2],
        "Microglia" = IBM_cols[3]
      )
    ))

pca_result <-  prcomp(t(UKBBN2_betas),        
                      retx = TRUE,
                      center =TRUE,
                      scale = FALSE,
                      rank. = 5)
scores_only_all <- as.data.frame(pca_result$x)

#PCA results are combined with phenotype data
pca_result_sample_info <- scores_only_all %>%
  mutate(ID = row.names(scores_only_all)) %>%
  inner_join(UKBBN2_pheno, by = c("ID" = "Basename"))

(PCA_celltype_UKBBN2 <- ggplot(pca_result_sample_info, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = `cell type:ch1` ), alpha = 0.7) +
    coord_fixed(1) +
    labs(
      x = paste0("PC1 (",round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 2), "%)"),
      color = "Cell Type") + 
    scale_color_manual(
      values = c(
        "Oligodendrocytes" = IBM_cols[1],
        "Astrocyte-enriched" = IBM_cols[4],
        "Neurons" = IBM_cols[2],
        "Microglia" = IBM_cols[3]
      )
    ))

#### combine all PCAs into one supplementary figure
# extract the legend
legend_grob <- cowplot::get_legend(PCA_celltype_BDR + theme(legend.position = "right"))

# combine the plots and use the empty space for the legend
final <- cowplot::plot_grid(
  PCA_celltype_UKBBN2 + theme(legend.position = "none"),
  PCA_celltype_UKBBN2_NC + theme(legend.position = "none"),
  PCA_celltype_BDR + theme(legend.position = "none"),
  legend_grob,
  ncol = 2,
  labels = c("A", "B", "C", ""),
  rel_widths = c(1, 1)
)

# save the figure
ggsave("Figures/S1.svg", device = "svg", width = 160, height = 160, units = "mm")
