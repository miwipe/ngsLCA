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

R script ngsLCA_interpret.R processes and visualizes outputs from ngsLCA.

It will detecte and install required packages automaticlly before running, and was developed and tested under R version 3.6.1.

Bugs and/or suggestions to wyc661217@gmail.com

## Input files

For inputing, this script requires all lca files (*.lca) copied into a working directory.

Two optional inputs:

1) lca files for negetive laboratory control samples can be copied into a sub-directory named "blanks" under the working directory. If provided, the taxa detected in control samples will be removed from real samples.

2) A csv metadata file matches the file names to sample names, ages, locations, etc. if provided, the supplied metadata will be illustrated in the results instead of file names. The metadata file format should be a tab separated two columns flat text, with first column covering lca file names, and second column supplying the metadata that will be illustrated. An example "metadata.txt" can be found under the R folder.

## Parameters

An example for running ngsLCA_interpret.R:

Rscript path_to_script/ngsLCA_interpret.R path="working_directory/" func = c("NMDS", "group", "rarefy", "heatmap") thr1=2 thr2=5 metadata="path_to_metadata/metadata.txt" taxa.re = c(1，131567) group.name = c("2:Bacteria", "33090:Viridiplantae") top.abundance = 30

Parameters:

path -- working directory containing all lca files

func -- functions that will be performed; default: NMDS, group, rarefy, heatmap; other option: stratplot (recommend            when metadata are ages) 

thr1 -- minimum reads number required for a taxon in each sample, integer, default: 2

thr2 -- minimum summed reads number required across all samples for a taxa, default: 5
      
metadata -- full path to your metadata, optional

taxa.re -- a list of NCBI taxaID representing the taxa that will be removed from final results, taxaID can be found at https://www.ncbi.nlm.nih.gov/Taxonomy

sample.re -- a list of lca file names that will not be included in the final results, e.g. sample.re = c("file1.lca","file5.lca")

group.name -- higher taxonomic ranks that will be used for grouping taxa, format: "NCBI taxaID:Scientific name"; default: "2:Bacteria", "33630:Alveolata", "33682:Euglenozoa", "4751:Fungi", "33208:Metazoa", "33090:Viridiplantae", "10239:Viruses"

thr3 -- minimum percentage of the reads of a taxon to the total reads numbers of the group, range from 0 to 1, default: 0

top.abundance -- how many most abundant taxa will be illustrated in figs, default 100















