# Cell-type Enrichment Analysis for Methylation (CEAM)
The developed method for cell type enrichment using methylation data is an over-representation analysis approach. This repository contains the scripts used for constructing the cell type-specific CpG sets, testing the method with simulated data and applying it to various EWAS results. The data used for establishing this method are not publicly avaialable, however the resulting CpG sets and enrichment results are available within this repository or in the supplementary materials found here (bioXriv link). The R scripts used for this method can be found in the Scripts folder.

The current version of CEAM was constructed on the on the four main brain cell types: neurons, microglia, oligodendrocytes, and astrocytes, using nuclei-sorted samples of neurologically elderly brains. This makes the current version of the CEAM highly specific for use cases in elderly brains such as neurodegenerative disorders. However, phenotypes far outside the scope of elderly brains and cancer may be not applicable to the current function of CEAM. Moreover, the brain samples used in the construction of this method were obtained from the prefrontal cortex and cingulate gyrus. Therefore, applying the tool to results of other brain regions may yield unexpected or spurious results. Although our application to EWAS results have shown enrichment results to be plausible within the light of current disease knowledge, despite these EWAS results being obtained from different brain regions.

# Citation
If you use the CEAM tool or the cell type-specific CpG sets please reference:

Müller J, Laroche V, Imm J, Weymouth L, Harvey J, Smith AR, et al. A Cell Type Enrichment Analysis Tool for Brain DNA Methylation Data (CEAM) [Internet]. bioRxiv; p. 2025.07.08.663671. https://doi.org/10.1101/2025.07.08.663671

# Quick Start
For a fast and easy cell type enrichment analysis we provide an app hosted here: (https://um-dementia-systems-biology.shinyapps.io/CEAM/)
Or download the Enrichment_Analysis_App.zip folder to run it locally

Steps:
1. Paste or import your list of differentially methylated probe identifiers 
2. Paste or import your background from which probe identifiers you have obtained these DMPs
3. Select the specificity level of CpG set you want to use (further explained below)
4. Extract the relevant resulting figures or table

# Workflow of the Tool
CEAM is based on an over-representation analysis framework. Therefore, cell type-specific CpG sets were constructed focusing on intermediately methylated CpGs within each cell type. Moreover, multiple sets per cell type were constructed as different specificity levels:

High-specificity: intermediate methylation value in only one cell type

Medium-specificity: intermediate methylation value in multiple cell types

Low-specificity: methylation value significantly different from all other cell types

These CpG sets were constructed with increasing coverage but decreasing cell type-specificity. The tool performs a hypergeometric test comparing the user's input CpGs to the user-specified CpG set of the desired specificity level and returns the resulting odds ratio, P-value and q-value (adjusted p-value for multiple testing).

# Interpretation of the Results
Since the CpG sets for each cell type are constructed with a focus on intermediately methylated cytosines, the interpretation remains straightforward. Enrichment in a cell type indicates that from the differentially methylated CpGs more than by chance are overlapping with the intermediately methylated CpGs within a cell type. Since intermediately methylated CpGs have been shown to have a higher cell-to-cell and inter-individual variability, it is likely that the observed signal in your bulk data is driven by changes in these intermediately methylated CpGs within the enriched cell type(s).
To summarize: enrichment of a cell type indicates that observed bulk changes are likely occurring within that cell type.


Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg