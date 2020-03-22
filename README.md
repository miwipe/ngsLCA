# ngsLCA
[![Build Status](https://travis-ci.org/miwipe/ngsLCA.svg?branch=master)](https://travis-ci.org/miwipe/ngsLCA)

This is the official development repository for ngsLCA (next generation sequence Least Common Ancestor algorithm). This package includes two modules, of which the first provides a fast and flexible taxonomic classification of DNA reads aligned to reference databases containing multiple organisms. The classification builds upon the NCBI taxonomy and performs a naïve least common ancestor assignment for reads with multiple alignments against different references. It is a commandline based tool that outputs a text file in "lca" format. 

An Rscript as the second module is provided for a quick transformation of the "lca" files to tables in different formats after filtering, e.g. a regular tab seperated table, krona and megan compatible formats. The output tables are also split into different kingdoms (or user-defined taxonomic groups) and taxonomic levels. In addition, it outputs 
heatmaps, barplots and stratplots as well as NMDS and rarefaction analysis for a quick overview of the dataset. If laboratory controls have been sequenced and are provided it can also subtract contamination taxa from the output. 

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
For a quick test of the ngsLCA, alignment files in bam format as inout can be found in the folder: "bam_files".

To generate alignment files from your own data (in bam/sam) please follow this quick guide on how to prepare your own data:

1. Download raw sequencing data, example fastq-files can be found in the fastq folder (It is assumed that all fastq-files have been demultiplexed, trimmed and quality controlled). 

2. Download a database of your own choice. The toolit is built upon the NCBI taxonomy therefore requires the reference database(s) complying to the NCBI format, i.e. fasta header should contain accession ID as the first string that appears in the NCBI access2taxID file, and the corresponded taxID is present in NCBI taxonomy dmp files. This could be the NCBI-nt and NCBI-RefSeq, for example NCBI-Refseq plastid database: 

```
mkdir refseq_plastids;
cd refseq_plastids;
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz;
gzip -d *;
cat *.fna > plastids.fa;
rm *.fna;
bowtie2-build --threads 5 plastids.fa plastids 
```

3. Align your trimmed and quality checked reads against the database 
```
bowtie2 --threads 10 -k 5000 -x refseq_plastids/plastids -U fastq_file.fq --no-unal | samtools view -bS - > file_name.database_name.bam
```

4. If more than one database have been used as reference all resulting bam files needs to be merged and sorted using samtools (important as the LCA assumes that all unique readIDs are aligned next to each other). See example below:

```
samtools merge -@ 10 -n merged.out.bam *.bam
samtools sort -n -T /TMP_folder/ -O bam -o file.sort.bam -@ 5 tmp.bam.merged -m 2G
```

# Running the ngsLCA first module
## Downloading resource files for program from NCBI 
```
mkdir ncbi_tax_dmp;
cd ncbi_tax_dmp/;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip;
unzip new_taxdump.zip;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz;
gunzip nucl_gb.accession2taxid.gz;
```

## Running the main ngsLCA program
The first module takes into account a chosen similarity interval between each read and its reference in the generated bam/sam file. It has been tested downloading the nt, refseq as well as individual genomes and fasta sequences. The similarity can be set as either an edit distance [-editdist[min/max]] eg. number of mismatches between the read and each alignment to a reference genomes or as a similarity distance [-simscore[low/high]] eg. a percentage of mismatches between the read and each alignment to a reference genomes.

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
The similarity score function works exactly the same way as for the edit distance, with the exception that the similarity to the reference is now calculated as a percentage and hence takes the length of the read into account. Example for running the ngsLCA algorithm with similarity scores (95-100% similarity) [-simscore[low/high]] 0.95-1.0. 
```
ngsLCA/ngsLCA -simscorelow 0.95 -simscorehigh 1.0 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile.ss095to1
```

## The .lca output file format
The resulting file (.lca) from the first module is a regular text file in which each line contains a unique read that aligned to a reference with the given simlarity speficified when executing ngsLCA. 

The file is tab seperated and the first column contains the read metadata which is subdivided in a colon seperated style and the first is the read ID (colon seperated column 1-7), the query sequence (csc 8), the lenght of the sequence (csc 9), the number hits to a reference in the database(s) (scs 10). 

The second column contain the lowest taxonomic node the read has been assigned to, seperated by colon is the 'NCBI taxID':'taxon name':'taxonomic level' assigned to. Following columns contain the taxonomic path higher in the NCBI taxonomy for each assignment.  


# Running the ngsLCA second module

The script ngsLCA_interpret.R is supplied in the folder "R_script", which filters, splits and visualizes outputs of the first module. It assumes each input lca file representing a sample and will generate a folder in the working directory containing all results. 


## Input files

This script requires all .lca files to be located in the same working directory. 

An example for running ngsLCA_interpret.R:
```
Rscript ngsLCA/R_script/ngsLCA_interpret.R path="working_directory" metadata="path_to_metadata/metadata.txt"
```

## Parameters

Run the script without any input to show the default value for each parameter:
```
Rscript ngsLCA/R_script/ngsLCA_interpret.R
```

Parameters:

path -- working directory containing the lca files that will be processed
###### This is the only compulsory parameter for the script

path_blank -- working directory containing the lca files for all blank controls
###### If supplied, a contamination list will be generated by merging all taxa deteced in these control samples, and the listed taxa will be substracted (optional) from the uotput results.
  
output -- name of the output folder under the working directory (path)

task -- functions that will be performed, the "pre-process" and "filter" need to be first performed before running other    functions.
  pre-process -- count the reads number for each detected taxon to form a taxon-count matrix for each input lca file; merging all taxon-count matrixes to generate a combined taxonomic profile.
  filter -- filter the taxonomic profile by user-defined thresholds.
  de-contamination -- substract the taxa in the contamination list from the taxonomic profile, if path_blank supplied.
  group -- split the taxonomic profile into different kingdoms (or user-defined taxonomic groups).
  rank -- split the taxonomic profile into different taxonomic ranks.
  count -- count the reads number and taxa number after each filter.
  megan -- generate file for MEGAN input.
  krona -- generate file for krona input.
  heatmap -- generate heatmaps.
  barplot -- generate barplots.
  stratplot -- generate stratplots
  rarefy -- perform the random rarefaction analysis.
  NMDS -- perform the NMDS rarefaction analysis.

metadata -- path to your metadata
  The supplied metadata should be in a tab (\t) separated format, and containing three columns. In which the first column contains a list of all lca file names, second column should contain the desired naming of the samples in the illustrations, and third column should be a numeric value that can order samples (age, depth or rank, interval are not allowed). An example "metadata.txt" can be found under the "R_script" folder.

threshold.1 -- minimum reads number required for confirming a taxon in each sample

threshold.2 -- minimum reads percentage (to the total reads number of the sample) required for confirming a taxon in each sample, ranging from 0 to 1

threshold.3 -- minimum summed reads number across all samples required for confirming a taxon in the combined taxonomic profile

threshold.1_blank -- same as threshold.1, but for filtering the blank controls; will using the value of threshold.1 if not specified

threshold.2_blank -- same as threshold.2, but for filtering the blank controls; will using the value of threshold.2 if not specified

threshold.3_blank -- same as threshold.3, but for filtering the blank controls; will using the value of threshold.3 if not specified
 
remove.taxa -- a list of NCBI taxaID indicating the taxa that will be removed from final results

remove.sample -- a list of file names indicating the samples that will be removed from the final results

group.name -- higher taxonomic units that will be used for grouping the taxa, need to be as the scientific names in NCBI taxonomy

rank.name -- taxonomic ranks that will be used for clustering the taxa

threshold.perGroup -- minimum reads percentage (to the total reads number of each group) required for confirming a taxon in each sample, ranging from 0 to 1

top.abundance -- number of most abundant taxa that will be illustrated in figs

NMDS_trymax -- maximum numbers of random starts in search of convergent solutions for NMDS


## Results

intermediate -- intermediate files, containing the taxonomic profile after each filtering
taxonomic_profiles -- tab separated tables for all the cleaned taxonomic profiles
counts -- reads and taxa number statistics for each sample
megan -- generated files for MEGAN input
krona -- generated files for krona input
heatmap -- generated heatmaps
barplot -- generated barplots
stratplot -- generated stratplots
rarefaction -- random rarefaction curves
NMDS -- NMDS analysis output data and figures
