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

The R script ngsLCA_interpret.R processes and visualizes outputs from ngsLCA. It will generate a directory "R_results" in your appointed working directory. The script will detecte and install required packages automaticlly before running. 

Developed and tested under R version 3.6.1. Bugs and/or suggestions to wyc661217@gmail.com

## Input files

For inputing, the script requires all lca files (*.lca) copied into a working directory.

Two optional inputs:

1) lca files for negetive laboratory control samples can be copied into a sub-directory named "blanks" under the working directory. If provided, the taxa detected in control samples will be removed from real samples.

2) A csv metadata file matches the file names to the metadata and orders for illustrating your samples. if provided, the supplied metadata will be shown in the results instead of file names. The metadata file format should be a tab (\t) separated three columns flat text, with first column covering lca file names, second column supplying the metadata that will be illustrated, and third column for the order of files. An example "metadata.txt" can be found under the R folder.

## Parameters

An example for running ngsLCA_interpret.R:

Rscript path_to_script/ngsLCA_interpret.R path="working_directory/" metadata="path_to_metadata/metadata.txt"

Parameters:

path -- working directory containing all lca files

func -- functions that will be performed; default: func="NMDS, group, rarefy, heatmap"; other option: stratplot (recommend            when metadata are ages) 

thr1 -- minimum reads number required for a taxon in each sample, integer, default: thr1=2

thr2 -- minimum summed reads number required across all samples for a taxa, default: thr2=5
      
metadata -- full path to your metadata, optional

taxa.re -- a list of NCBI taxaID representing the taxa that will be removed from final results, taxaID can be found at https://www.ncbi.nlm.nih.gov/Taxonomy; taxa="1,131567" will romove "root" with taxaID 1 and "cellular organisms" with taxaID 131567

sample.re -- a list of lca file names that will not be included in the final results, sample.re = c("file1.lca","file5.lca")

group.name -- higher taxonomic ranks that will be used for grouping taxa, format: "NCBI taxaID:Scientific name"; default: group.name="10239:Viruses,2157:Archaea,2:Bacteria,4751:Fungi,33090:Viridiplantae,33317:Protostomia,33208:Metazoa"

thr3 -- minimum percentage of the reads for a taxon to the total reads numbers of the group, range from 0 to 1, default: thr3=0

top.abundance -- how many most abundant taxa will be illustrated in figs, default top.abundance=50















