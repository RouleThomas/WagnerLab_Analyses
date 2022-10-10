#!/bin/R

# run this script using the following command:
# Rscript [path/to/script]/install_R_packages.R

# install packages from CRAN
install.packages("argparse",repos = "http://cran.us.r-project.org")
install.packages("tidyverse",repos = "http://cran.us.r-project.org")
install.packages("ggplot2",repos = "http://cran.us.r-project.org")

# install packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("DESeq2")

#shrinkage methods
BiocManager::install("apeglm")
BiocManager::install("ashr")

# Print out package info
print("PACKAGE VERSIONS:")
print(paste("argparse"))
print(packageVersion("argparse"))
print("tidyverse")
packageVersion("tidyverse")
print("ggplot2")
packageVersion("ggplot2")
print("DESeq2")
packageVersion('DESeq2')
print('apeglm')
packageVersion('apeglm')
print('ashr')
packageVersion('ashr')
