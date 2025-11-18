source("Scripts/Functions/Functions.R")

load("Data/UKBBN2_betas.RData")
load("Data/BDR_betas.RData")

#### EPIC background
# intersect both cohorts after preprocessing
EPIC_bkgd <- (intersect(rownames(BDR_betas), rownames(UKBBN2_betas)))
write.table(EPIC_bkgd, file = "CpG sets/EPIC_background.txt", sep = "\n", row.names = F, col.names = F)

#### 450K background

# load 450K probes
load("Data/background_450K.RData") # all probes on a 450K array

# filter out problematic probes such as corss hybridising and SNPs
crosshyb <- read.csv("Data/CrossHybridisingProbesPriceORWeksberg.csv")
background_450 <- background_450[!background_450 %in% crosshyb$x]

### Removing SNP CpG sites
SNP1 <- read.csv("Data/snp_chyb_Af_out.csv")
SNP2 <- read.csv("Data/snp_chyb_EU_out.csv")
SNPall <- union(SNP1$Probe, SNP2$Probe)

background_450 <- background_450[!background_450 %in% SNPall]

#### Removing X and Y chromosomes probes
manifest <- read.csv("Data/humanmethylation450_15017482_v1-2.csv", header = T, skip = 7)
background_450 <- background_450[!background_450 %in% (manifest[manifest$CHR %in% c("X", "Y"), "IlmnID"])]

# intersect the EPIC background with every probe remaining in 450K
background_450 <- intersect(background_450, EPIC_bkgd)
write.table(background_450, file = "CpG sets/background_450k.txt", sep = "\n", row.names = F, col.names = F)


