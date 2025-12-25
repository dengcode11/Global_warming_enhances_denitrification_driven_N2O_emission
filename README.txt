Dengm 2025-Global warming enhances denitrification-driven N2O emission from plastisphere and epiphyton

1. Scripts for bioinformatics analysis
2. R scripts for statistical analysis and plotting figures
3. Data for plotting Figures 1-3 and Figure 5. All the figures are created by Graphpad Prism 9.0.0, R software (v4.2.0) and Adobe illustrator 2022
4. Climate data of five widely-used Earth System Models (ESM) under the SSP 3-RCP 7.0 scenario
5. Sequences for phylogenetic trees of nosZ


# 1.1 Projected increase in denitrification rates and denitrification-driven N2O emission potential

## Installation

R        (v4.2.0) (https://www.r-project.org/) 
R-studio (v2023.09.0+463)   (https://posit.co/products/open-source/rstudio/) # IDE for R

### Required packages
ncdf4     (v1.22)    # For reading NetCDF files
raster    (v3.6-26)  # For handling raster (gridded) data
tidyverse (v2.0.0)   # For data manipulation and visualization
lubridate (v1.9.3)   # For date-time operations
patchwork (v1.3.2)   # For combining multiple plots
rstatix   (v0.7.2)   # For statistical tests
ggpubr    (v0.6.0)   # For publication-ready plots
ggsignif  (v0.6.4)   # For adding statistical significance annotations to plots

### Usage
The climate data are provided in fold named "climate data of five widely-used ESM under the SSP 3-RCP 7.0 scenario". This data were download from https://aims2.llnl.gov/search/cmip6, which is the official data portal of the Coupled Model Intercomparison Project Phase 6 (CMIP6), hosting the core climate model simulations used in the latest IPCC assessment reports.

For detailed guidance and functional comments on the code, please refer to the file named "End-Century_N2O_Emission_Increases_A_Q10-Based_CMIP6_Multi-Model_Analysis.R", where we have added detailed comments to each line of code.


### Main procedures for analysis
1. Define Key Parameters       # Q10 value
2. Define File Paths           # File Paths of climate data
3. Define Helper Functions
4. Process All Model Data
5. Calculate ΔR% for Each Q10 Value
6. Prepare Boxplot Data (2090-2100)
7. Perform Statistical Analysis
8. Create Boxplot
9. Calculate Ensemble Statistics (for main trend plot)
10. Create Main Trend Plot
11. Save All Data and Results


# 1.2 Statistical analysis for denitrification rate, denitrification-driven N2O, and Ea
## Installation

R        (v4.2.0) (https://www.r-project.org/) 
R-studio (v2023.09.0+463)   (https://posit.co/products/open-source/rstudio/) # IDE for R
### Required packages
tidyr    (v1.3.0)
dplyr    (v1.1.4)
rstatix  (v0.7.2)
ggplot2  (v3.5.2)
ggpubr   (v0.6.0)
effsize  (v0.8.1)
coin     (v1.4-3)
car      (v3.1-3)


### Usage
For detailed guidance and functional comments on the code, please refer to the file named "Statistical_analysis_for_denitrification_N2O_Ea.R", where we have added detailed comments to each line of code.

### Main procedures for analysis
1. DATA INPUT    
2. INDEPENDENT SAMPLES COMPARISON (Differences between groups for the same variable)           
3. PAIRED SAMPLES COMPARISON (Differences between variables within the same group)
4. RESULTS SUMMARY
5. MULTIVARIATE VISUALIZATION
6. SAVE RESULTS



# 1.3 Metagenomic Analysis Project

## Required Packages
For the installation method, please refer to the following URL.

fastp (v0.23.4) (https://github.com/OpenGene/fastp)

MEGAHIT (v1.2.9) (https://github.com/voutcn/megahit)

Bowtie2 (v2.5.1) (https://github.com/BenLangmead/bowtie2)

samtools (v1.19) (https://github.com/samtools/samtools)

MetaBAT2 (v2.15) (https://bitbucket.org/berkeleylab/metabat)

dRep (v3.4.0) (https://github.com/MrOlm/drep)

CheckM (v1.2.2) (https://github.com/Ecogenomics/CheckM)

GTDB-Tk (v2.4.0) (https://github.com/Ecogenomics/GTDBTk)

CoverM (v0.6.1) (https://github.com/wwood/CoverM)

Prodigal (v2.6.3) (https://github.com/hyattpd/Prodigal)

HMMER (v3.4) (https://github.com/EddyRivasLab/hmmer)

bbmap (v39.01) (https://github.com/BioInfoTools/BBMap)

Muscle (v5.1) (https://github.com/rcedgar/muscle)

trimAl (v1.4.1) (https://github.com/inab/trimal)

IQ-TREE (v2.3.1) (https://github.com/iqtree/iqtree2)

METABOLIC (v4.0) (https://github.com/AnantharamanLab/METABOLIC) 


## Usage
After installing the required packages, you can run the metagenomic analysis code as described in the file named "Metagenomic_analysis.txt". For detailed guidance and functional comments on the code, please refer to the file named "Metagenomic_analysis.txt".

### Main procedures for analysis
 1. Quality control
 2. Metagenomic Assembly
 3. Aligning sequencing data to a reference genome
 4. Processing SAM/BAM files
 5. Metagenomic binning
 6. Metagenomic dereplication
 7. Evaluating the quality of metagenomic bins
 8. Metagenomic classification
 9. Calculating genome coverage
 10. Metagenomic functional annotation
 11. Performing gene prediction on contigs
 12. Performing functional search on the predicted protein sequences
 13. Extracting sequences with specific IDs from a nucleic acid sequence file
 14. Build an index for the extracted gene sequence file
 15. Calculate the RPKM values of the alignment results
 16. Phylogenetic tree analysis

#### Notes
 Ensure that all input files (e.g., Rawreads_1.fq.gz, Rawreads_2.fq.gz, final.contigs.fa, Gene.hmm) are correctly placed in the specified directories. 

 Adjust the paths in the commands as needed based on your project structure. 

