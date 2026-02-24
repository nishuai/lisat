[![Build Status](https://travis-ci.org/lima1/PureCN.svg?branch=master)](https://travis-ci.org/lima1/PureCN)
[![Coverage](https://img.shields.io/codecov/c/github/lima1/PureCN.svg)](https://codecov.io/gh/lima1/PureCN)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0) 

# 💡 Lisat

<!-- badges: start -->
<!-- badges: end -->

## Overview

**`lisat`** is a comprehensive R toolkit designed for the analysis of longitudinal virus integration site data. It streamlines the entire workflow of integration site analysis, from data cleaning and quality control to statistical modeling and rich visualization. With support for simple input formats, `lisat` provides a user-friendly and powerful suite of functions for researchers investigating viral integration sites, clonal tracking, and gene therapy safety.

## Key Features

- **Genomic Feature Annotation**: Automatically annotate integration sites with genomic features. Easily check for overlaps with critical genomic regions including enhancers, promoters, safe harbors, adverse event (AE) genes, cancer genes, and immune-related genes.
- **Integration Site Analysis (CIS)**: Identify Common Integration Sites (CIS) to discover regions with recurrent integrations and analyze overall chromosome distributions.
- **Longitudinal & PMD Analysis**: Track clonal dynamics over multiple timepoints using Population Matching Distribution (PMD) analysis. Evaluate clonal richness and evenness across different patient samples.
- **Clonal Dominance Analysis**: Identify and analyze potential dominant clones using cumulative distribution fitting models.
- **Rich Visualization**: Built-in plotting functions for creating high-quality, publication-ready visualizations, including treemaps, cumulative curves, region count pie charts, and chromosome ideograms.

## Installation

You can install the development version of `lisat` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("nishuai/lisat")
```


### Dependencies

For full annotation capabilities, ensure the following Bioconductor packages are installed:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
```

## Quick Start

Here is a basic example showing how to validate your raw data and perform an initial analysis:

``` r
library(lisat)

# 1. Prepare your raw integration site data
# (Requires columns: Sample, SCount, Chr, Locus)
head(IS_raw)

# 2. Validate the data structure
check_validity <- validate_IS_raw(IS_raw)

# 3. Annotate Genomic Features
# Requires TxDb.Hsapiens.UCSC.hg38.knownGene and org.Hs.eg.db
IS_annotated <- get_feature(IS_raw)
IS_annotated <- Enhancer_check(IS_annotated)
IS_annotated <- Promotor_check(IS_annotated)
IS_annotated <- Safeharbor_check(IS_annotated)

# 4. Identify Common Integration Sites (CIS)
CIS_top <- CIS(IS_raw = IS_annotated, connect_distance = 50000)
CIS_overlap(CIS_data = CIS_top, IS_raw = IS_annotated)

# 5. Longitudinal Analysis
# Requires a Patient_timepoint metadata dataframe
PMD_data <- pmd_analysis(IS_raw = IS_annotated, Patient_timepoint = Patient_timepoint)
plot_richness_evenness(PMD_data = PMD_data)
```

For a comprehensive guide, please refer to the package vignette:
``` r
vignette("lisat-intro", package = "lisat")
```

## Citation

If you use `lisat` in your research, please cite our preprint:

> Ni, S. et al. (2025). *[Insert exact title from bioRxiv here]* bioRxiv. DOI: [10.64898/2025.12.20.695672v1](https://www.biorxiv.org/content/10.64898/2025.12.20.695672v1)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
