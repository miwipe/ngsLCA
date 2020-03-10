# ngsLCA
[![Build Status](https://travis-ci.org/miwipe/ngsLCA.svg?branch=master)](https://travis-ci.org/miwipe/ngsLCA)

This is the official development repository for ngsLCA (next generation sequence Least Common Ancestor algorithm). This package provides a fast and flexible taxonomic classification of DNA reads aligned to a reference database containing multiple organisms. The classification builds upon the NCBI taxonomy and performs a naïve least common ancestor assignment for reads with multiple alignment against different references. It is a commandline based tool that outputs a text file whcih easily can be parsed and handled in R or like, for further statistical analysis. 

An Rscript is provided for a quick transformation of the lca output to tables which are provided in different formats e.g. a regular comma seperated table as well as krona and megan compatible input formats. The output tables is also split into different kingdoms and taxonomic levels. In addition, it output NMDS plots, rarefaction analysis, heatmaps and stratplots for a quick overview of the dataset. Lastly, if laboratory controls have been sequenced and are provided it can subtract contamination taxa from the final output. 

# Building ngsLCA
ngsLCA requires [HTSlib](https://github.com/samtools/htslib) which is a common library used for handling high-throughput sequencing data. You can install it as shown below or link to a previously-installed HTSlib when running make on ngsLCA.  


```
git clone https://github.com/SAMtools/htslib
git clone https://github.com/miwipe/ngsLCA
cd htslib
make
cd ../ngsLCA
make HTSSRC=../htslib
```
# Test dataset
For a quick test of the installation, alignment files in bam formats can be found in the folder: "bam_files"

To generate alignment files from your own data (in bam/sam) please follow this quick guide on how to prepare your own data:
1. Download raw sequencing data, (example fastq files can be found here "LINK TO FASTQ" (It is assumed that all fastq files have been demultiplexed, trimmed and quality controlled). 
2. Download a database of your own choice but based on the ncbi taxonomy (it is required that the fasta header in the database contains ncbi accession no. provided by ncbi and that this appears in the first field of the header). This could be the RefSeq plastid database: 

```
mkdir refseq_plastids;
cd refseq_plastids;
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz;
gzip -d *;
cat *.fna > plastids.fa;
rm *.fna;
bowtie2-build --threads 5 plastids.fa plastids 
```

3. Next, align your trimmed and quality checked reads against the database 
```
bowtie2 --threads 40 -k 5000 -x refseq_plastids/plastids -U fastq_file.fq --no-unal | samtools view -bS - > file_name.database_name.bam
```

If more than one database have been used as reference all resulting bam files needs to be merged and sorted using samtools (important as the LCA assumes that all unique readIDs are aligned next to each other). See example below:

```
samtools merge -@ 10 -n merged.out.bam *.bam
samtools sort -n -T /TMP_folder/ -O bam -o file.sort.bam -@ 5 tmp.bam.merged -m 5G
```

# Running Main c/c++ program
## Downloading resource files for program from NCBI 
```
mkdir ncbi_tax_dmp;
cd ncbi_tax_dmp/;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip;
unzip new_taxdump.zip;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz;
gunzip nucl_gb.accession2taxid.gz;
```

## Running ngsLCA 
The LCA algorithm is a naïve LCA that calculates the least common ancestor using the NCBI taxonomy for a multiple alignment file in bam format sorted by read ID. It takes into account a chosen similarity interval between each read and its alignments. The headers of the database therefore needs to follow the nomenclature from NCBI. It has been tested downloading the nt, refseq as well as individual genomes and fasta sequences. The similarity can be set as either an edit distance [-editdist[min/max]] eg. number of mismatches between the read and each alignment to a reference genomes or as a similarity distance [-simscore[low/high]] eg. a percentage of mismatches between the read and each alignment to a reference genomes.

Edit distance can be a number between 0-10, while the similarity score is a number between 0-1. 

Example for running the ngsLCA algorithm with 0 edit distances to reference.
```
ngsLCA/ngsLCA -editdistmin 0 -editdistmax 0 -names ncbi_tax_dmp/names.dmp.gz -nodes ncbi_tax_dmp/nodes.dmp.gz -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid_24april.gz -bam sorted_bam_file.bam -outnames outfile.ed0
```

Example for running the ngsLCA algorithm with an interval of edit distances (min 0, max 5) to reference. In this case the algorithm will perform an LCA on all alignment of a given read that have between 0 and 5 mismatches. NOTE! It will not discriminate between better or worse alignments within the defined interval.
```
ngsLCA/ngsLCA -editdistmin 0 -editdistmax 5 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile
```
If the editditmin is for example set to 2 the algorithm will only parse reads with the best alignment of 2 or more mismatches and not parse reads that have better alignments. This allows for a stepwise examination of the taxonomic profiles generated by allowing gradually more mismatches. 
```
ngsLCA/ngsLCA -editdistmin 2 -editdistmax 2 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile.ed2
```
The similarity score function works exactly the same way as for the edit distance, with the exception that the similarity to the reference is now calculated as a percentage and hence takes the length of the read into account. Example for running the ngsLCA algorithm with similarity scores (90-100% similarity) [-simscore[low/high]] 0.9-1.0. 
```
ngsLCA/ngsLCA -simscorelow 0.90 -simscorehigh 1.0 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile.ss09to1
```

## The .lca output file format
The resulting file from ngsLCA (.lca) is regular text file in which each line contains a unique read that aligned to a reference with the given simlarity speficified when executing ngsLCA. The file is colon seperated and the first column contains the read metadata which is subdivided in a colon seperated style and the first is the read ID (colon seperated column 1-7), the query sequence (csc 8), the lenght of the sequence (csc 9), the number hits to a reference in the database(s) (scs 10). 

The second column (tab seperated) contains the lowest taxonomic node the read has been assigned to, seperated by colon is the 'NCBI taxID':'taxon name':'taxonomic level' assigned to. Following columns contains the taxonomic path higher in the NCBI taxonomy for each assignment.  


# Visualizing results with R

The R script ngsLCA_interpret.R merges, filters, splits and visualizes outputs from ngsLCA (.lca files). It will generate a directory "R_results" in the working directory. The script will detect and install required packages automatically before running. 

This script was developed and tested under R version 3.6.1. Bugs and/or suggestions email wyc661217@gmail.com and mwpedersen@sund.ku.dk

## Input files

This script requires all .lca files to be located in the same working directory. 

Before executing the script two options should be considered:

First, the script can remove taxa found in laboratory controls by making a folder in the working directory called 'blanks' and hereafter copying in the corresponding lca files from the labotory control. If provided, the taxa found in the control samples will be removed from the output from the 'true' samples.

Secondly, a metadata file can be supplied to rename and reorder samples (according to age, depth or a desired rank) in the evetual illustrations. This metadata should be in a tab (\t) separated file format, and containing three columns. In which the first column should contain a list of all lca file names (excluding the filenames in the blank folder), second column should contain the desired naming of the samples in the illustrations, and third column should be a numeric value that can order samples (age, depth or rank, interval are not allowed). An example "metadata.txt" can be found under the bam folder.

## Parameters

An example for running ngsLCA_interpret.R:
```
Rscript path_to_script/ngsLCA_interpret.R path="working_directory/" metadata="path_to_metadata/metadata.txt"
```
Parameters:

path -- working directory containing all lca files

func -- functions that will be performed; default: func="NMDS, group, rarefy, heatmap"; other option: stratplot (only works when metadata are supplied, best as depths or ages) 

thr1 -- minimum reads number required for a taxon in each sample, default: thr1=2

thr2 -- minimum summed reads number required across all samples for a taxa, default: thr2=5
 
thr3 -- minimum percentage of the reads for a taxon to the total reads numbers of the group, range from 0 to 1, default: thr3=0    

metadata -- full path to your metadata, optional

taxa.re -- a list of NCBI taxaID representing the taxa that will be removed from final results, taxaID can be found at https://www.ncbi.nlm.nih.gov/Taxonomy, for example inputting taxa="1,131567" will remove "root" with taxaID 1 and "cellular organisms" with taxaID 131567

sample.re -- a list of lca file names that will not be included in the final results, sample.re="file1.lca,file5.lca"

group.name -- higher taxonomic ranks that will be used for grouping taxa, format: "NCBI taxaID:Scientific name"; default: group.name="10239:Viruses,2157:Archaea,2:Bacteria,4751:Fungi,33090:Viridiplantae,33208:Metazoa"

top.abundance -- how many most abundant taxa will be illustrated in figs, default: top.abundance=50
