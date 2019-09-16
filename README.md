# ngsLCA
[![Build Status](https://travis-ci.org/miwipe/ngsLCA.svg?branch=master)](https://travis-ci.org/miwipe/ngsLCA)

# Installation
git clone https://github.com/SAMtools/htslib

git clone https://github.com/miwipe/ngsLCA

cd htslib;make;cd ../ngsLCA;make HTSSRC=../htslib

# Running Main c/c++ program
## Finding resource files for program from NCBI
## Running ngsLCA

# Visualizing results with R script

The R script ngsLCA_interpret.R processes and visualizes outputs from ngsLCA

Developed and tested in R version 3.6.1

The scripts relied on some R packages. It will detecte and install required R packages automaticlly before running.

Bugs and/or suggestions to ycwang@bio.ku.dk or wyc661217@gmail.com

## Input files

For inputing, this script requires all lca files (*.lca) copied into a working directory

Two optional inputs:

1) lca files for negetive laboratory control samples can be copied into a sub-directory named "blanks" under the working directory. If provided, the taxa detected in these control samples will be removed from real samples.

2) A metadata csv file matches the file names to sample names, ages, locations, etc. if provided, the supplied metadata will be illustrated in the results instead of file names. The metadata file accepte a comma separated two columns flat text, with first column covering all lca file names, and second column supplying the metadata will be illustrated. An example "metadata.txt" can be found under the R folder.

## Parameters

An example for running ngsLCA_interpret.R:

Rscript path_to_script/ngsLCA_interpret.R path="working_directory/" func = c("NMDS", "group", "rarefy", "heatmap") thr1=2 thr2=5 metadata="path_to_metadata/metadata.txt" taxa.re = c("1:root") group.name = c("2:Bacteria", "33090:Viridiplantae") top.abundance = 30

Parameters:

path -- working directory containing all lca files

func -- functions that will be performed; default: NMDS, group, rarefy, heatmap; other option: stratplot (recommend only            when metadata are ages) 

thr1 -- minimum reads number required for a taxon in each sample, integer, default: 2

thr2 -- minimum summed reads number representing a taxa across all samples, default: 5
      
metadata -- full path to your metadata, optional

taxa.re -- a list of taxa that will be removed from final results, format: "NCBI taxaID:Scientific name" e.g. "71240:eudicotyledons"

sample.re -- a list of samples that will be removed from the final results, file names with suffix ".lca" removed, e.g. sample.re = c("file1","file5")

group.name -- higher taxonomic ranks that will be used for grouping the taxa, format: "NCBI taxaID:Scientific name"; default: "2:Bacteria", "33630:Alveolata", "33682:Euglenozoa", "4751:Fungi", "33208:Metazoa", "33090:Viridiplantae", "10239:Viruses"

thr3 -- minimum percentage of the reads number of a taxon to the total reads numbers of the group, range from 0 to 1, default: 0

top.abundance -- how many most abundant taxa will be illustrated in figs, default 100















