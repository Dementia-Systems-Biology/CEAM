source("Scripts/Functions.R")
#### Required packages
check_and_load_libraries(c(
  "wateRmelon",
  "sva",
  "matrixStats",
  "ggplot2",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "RPMM",
  "dplyr",
  "tidyr"
))

#### Reading the data
# meta data
meta <- read.csv("UKBBN2_pheno.csv", header = T)

# reading idat files into extended rgSet object
rgSetEXT <- read.metharray.exp(base = "UKBBN2_FANS/", extended = T, targets = meta) # this function automatically recognizes the "basename" column

#### check bisulfite conversion #### 
bsc <- bscon(rgSetEXT) # all samples have good bsc

# transform data into dataframe for plotting
bsc_data <- data.frame(
  value = as.numeric(bsc),
  name = names(bsc)
)

# label the samples with value < 80
bsc_data$label <- ifelse(bsc_data$value < 80, bsc_data$name, NA) 

ggplot(bsc_data, aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "#01a08a", color = "black", alpha = 1) +
  geom_text(aes(label = label, y = 0), hjust = -0.5,angle = 90, color = "#b41921") +  # Add labels
  labs(title =  "Histogram of BS conversion", x = "bisulfite conversion percentage", y = "count")

# remove samples with bsc less than 80% (which were none)
rgSetEXT <- rgSetEXT[, !sampleNames(rgSetEXT) %in% bsc_data[bsc_data$value < 80, "name"]] 

#### Outlier detection ####
# detect outliers with IQR and PCOut algorithm
outliers <- outlyx(rgSetEXT, plot=FALSE)
# remove the samples that were deemed outliers by this method
rgSetEXT = rgSetEXT[, sampleNames(rgSetEXT) %in% rownames(outliers[outliers$outliers == FALSE,])] # 75 samples left
print(outliers[outliers$outliers == TRUE, ]) 


# Expression and outlier detection based on quality of the samples / CpG sites
rgSetEXT_pfilter <- pfilter(rgSetEXT)

# predicting sex
estimatesex = wateRmelon::estimateSex(getBeta(rgSetEXT_pfilter), do_plot=T) 
# merge estimated sex and observed sex information to compare
rownames(meta) <- meta$Basename
estimatesex_merged <- merge(estimatesex, meta, by = 'row.names', all.x = T, all.y = F)
# reformat Sex values
estimatesex_merged$Sex <- ifelse(estimatesex_merged$Sex == "M", "Male", "Female")

# find which do not match between predicted and observed and remove them
mismatch_sex <- estimatesex_merged[which(!estimatesex_merged$predicted_sex == estimatesex_merged$Sex), "Row.names"]

# remove the samples with mismatched sex
rgSetEXT_pfilter <- rgSetEXT_pfilter[, !sampleNames(rgSetEXT_pfilter) %in% mismatch_sex] 


#### Noob normalisation (dye-bias correction)
mSetEPIC_NOOB <- preprocessNoob(rgSetEXT_pfilter, offset = 15, dyeCorr = TRUE, verbose = FALSE,
                                dyeMethod="single")


#### BMIQ normalisation (intra-sample)
mset_betas_lm <- BMIQ(mSetEPIC_NOOB)

#### Removing all incomplete CpG sites
mset_betas_lm = mset_betas_lm[complete.cases(mset_betas_lm),]

#### visualization 
# pre-normalization betas
density_raw <- gather(as.data.frame(getBeta(rgSetEXT_pfilter)))
(density_raw_plt <- ggplot(density_raw) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Pre-normalization") +
  scale_color_viridis_d() +
  theme(legend.position =  "none"))

# post-normalization betas
density_nor <- gather(as.data.frame(mset_betas_lm))
(density_nor_plt <- ggplot(density_nor) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Post-normalization") +
  scale_color_viridis_d() +
  theme(legend.position =  "none"))

#### ComBat batch correction
# remove the meta data of the sample removed in the outlyx function
meta <- meta[meta$Basename %in% colnames(mset_betas_lm),]

if (all(meta$Basename == colnames(mset_betas_lm))){ # ensuring that samples are in the same order
  mset_betas_lm = ComBat(dat=mset_betas_lm, batch=as.character(meta$Institute), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  
}


### Removing SNP CpG sites
SNPall = read.table("SNPProbes_McCartney.txt", header = T) # obtained from McCartney DL, Walker RM, Morris SW, McIntosh AM, Porteous DJ, Evans KL. Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip. Genom Data. 2016 Sep;9:22–4. 
SNPallEURAF = SNPall[which(SNPall$EUR_AF >= 0.05 & SNPall$EUR_AF <= 0.95),]

#### Removing Cross-hybridising sites
crosshyb = read.table("CrossHydridisingProbes_McCartney.txt") # obtained from McCartney DL, Walker RM, Morris SW, McIntosh AM, Porteous DJ, Evans KL. Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip. Genom Data. 2016 Sep;9:22–4. 

#### Removing X and Y chromosomes probes
manifest_EPIC = read.csv("MethylationEPIC_v-1-0_B4.csv", header = T, skip = 7)
manifest_EPIC_XY = manifest_EPIC[manifest_EPIC$CHR %in% c("X", "Y"),]

to_remove_cpgs = unique(c(SNPallEURAF$IlmnID, manifest_EPIC_XY$IlmnID, crosshyb$V1))
#### From the 3 steps above you should remove 72240 CpGs in total

#### Removing Non-CpG sites
non_CpG_probes = rownames(mset_betas_lm[grep("ch.", rownames(mset_betas_lm)),])
mset_betas_lm = mset_betas_lm[!rownames(mset_betas_lm) %in% non_CpG_probes,]
dim(mset_betas_lm)

#### Removing Reference SNPs
ref_SNP_probes = rownames(mset_betas_lm[grep("rs", rownames(mset_betas_lm)),])
mset_betas_lm = mset_betas_lm[-(which(rownames(mset_betas_lm) %in% to_remove_cpgs)),]
mset_betas_lm = mset_betas_lm[!rownames(mset_betas_lm) %in% ref_SNP_probes,]

#### Removing Non-CpG sites
non_CpG_probes = rownames(mset_betas_lm[grep("ch.", rownames(mset_betas_lm)),])
mset_betas_lm = mset_betas_lm[!rownames(mset_betas_lm) %in% non_CpG_probes,]
dim(mset_betas_lm)

#### Performing PCA
pca_result <-  prcomp(t(mset_betas_lm),        
                      retx = TRUE,
                      center =TRUE,
                      scale = FALSE,
                      rank. = 5)
scores_only_all <- as.data.frame(pca_result$x)

#PCA results are combined with phenotype data
pca_result_sample_info <- scores_only_all %>%
  mutate(ID = row.names(scores_only_all)) %>%
  inner_join(meta, by = c("ID" = "Basename"))

(PCA_celltype <- ggplot(pca_result_sample_info, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = Cell_Type ), alpha = 0.7) +
    coord_fixed(1) +
    labs(
      x = paste0("PC1 (",round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
      color = "Cell Type"))

# rename the variables for later use in other scripts
UKBBN2_betas <- mset_betas_lm
UKBBN2_pheno <- meta # no meta data was removed because no samples were removed

save(UKBBN2_pheno, UKBBN2_betas, file = "UKBBN2_betas.RData")
