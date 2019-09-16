# ngsLCA
[![Build Status](https://travis-ci.org/miwipe/ngsLCA.svg?branch=master)](https://travis-ci.org/miwipe/ngsLCA)

# Installation
git clone https://github.com/SAMtools/htslib

git clone https://github.com/miwipe/ngsLCA

cd htslib;make;cd ../ngsLCA;make HTSSRC=../htslib

# Running Main c/c++ program
## Finding resource files for program from NCBI
## Running ngsLCA

# Processing and visualizing results using ngsLCA_interpret.R

This R script processes and visualizes outputs from ngsLCA 
Developed and tested in R version 3.6.1 under unix
Please report bugs and/or suggestions to ycwang@bio.ku.dk or wyc661217@gmail.com

### As inputing, this script requires all lca files (*.lca) copied into a working directory
# Two optional inputs:
# 1) lca files for negetive laboratory control samples can be  into a sub-directory named "blanks" under the working directory

# 2) and  .
# A metadata text (sample names, ages, locations, etc.) is optional as input. If provided, the reletive metadata will be illustrated in the results instead of file names
# see https://github.com/miwipe/ngsLCA for an example of metadata format
