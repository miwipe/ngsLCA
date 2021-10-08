# ngsLCA
[![Build Status](https://travis-ci.org/miwipe/ngsLCA.svg?branch=master)](https://travis-ci.org/miwipe/ngsLCA)

This is the official development repository for ngsLCA (next generation sequence Least Common Ancestor algorithm). This toolkit includes two modules, ngsLCA main program and ngsLCA R package. The ngsLCA main program provides a fast and flexible taxonomic classification of DNA reads aligned to reference databases containing multiple organisms. The classification builds upon the NCBI taxonomy and performs a naïve least common ancestor assignment for reads with multiple alignments against different references. It is a commandline based tool that outputs a text file in "lca" format. 

The ngsLCA R package provides functionality for processing "lca" files and illustrating the taxonomic profiles. It supplies R functions for quick transformation of the "lca" files to tables in different formats after filtering, e.g. a regular tab seperated table and megan compatible formats, for de-contamination if laboratory controls have been sequenced and are provided, and for spliting taxa into different kingdoms (or user-defined taxonomic groups) and taxonomic ranks. Functions are also provided for outputing heatmaps, barplots and stratplots as well as NMDS and rarefaction analysis for a quick overview of the dataset.

# Building ngsLCA
ngsLCA main program requires [HTSlib](https://github.com/samtools/htslib) which is a common library used for handling high-throughput sequencing data. You can install it as shown below or link to a previously-installed HTSlib when running make on ngsLCA.  


```
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/miwipe/ngsLCA
cd htslib
make
cd ../ngsLCA
make HTSSRC=../htslib
```

Released version of ngsLCA R package can be installed from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ngsLCA")
```

or from [github](https://github.com/wyc661217/ngsLCA) with:

``` r
remotes::install_github("wyc661217/ngsLCA")
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

4. If more than one database has been used as reference all resulting bam files needs to be merged and sorted using samtools (important as the LCA assumes that all unique readIDs are aligned next to each other). See example below:

```
samtools merge -@ 10 -n merged.out.bam *.bam
samtools sort -n -T /TMP_folder/ -O bam -o file.sort.bam -@ 5 tmp.bam.merged -m 2G
```

# Running the ngsLCA main program
## Downloading resource files from NCBI 
```
mkdir ncbi_tax_dmp;
cd ncbi_tax_dmp/;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip;
unzip new_taxdump.zip;
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz;
gunzip nucl_gb.accession2taxid.gz;
```

## Running the ngsLCA main program
The ngsLCA main program takes into account a chosen similarity interval between each read and its reference in the generated bam/sam file. It has been tested downloading the nt, refseq as well as individual genomes and fasta sequences. The similarity can be set as either an edit distance [-editdist[min/max]] eg. number of mismatches between the read and each alignment to a reference genomes or as a similarity distance [-simscore[low/high]] eg. a percentage of mismatches between the read and each alignment to a reference genome.

Edit distance can be a number between 0-10, while the similarity score is a number between 0-1. 

Example for running the ngsLCA main program with 0 edit distances to reference.
```
ngsLCA/ngsLCA -editdistmin 0 -editdistmax 0 -names ncbi_tax_dmp/names.dmp.gz -nodes ncbi_tax_dmp/nodes.dmp.gz -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid_24april.gz -bam sorted_bam_file.bam -outnames outfile.ed0
```

Example for running the ngsLCA main program with an interval of edit distances (min 0, max 5) to reference. In this case the algorithm will perform an LCA on all alignment of a given read that have between 0 and 5 mismatches. NOTE! It will not discriminate between better or worse alignments within the defined interval.
```
ngsLCA/ngsLCA -editdistmin 0 -editdistmax 5 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile
```
If the editditmin is for example set to 2 the algorithm will only parse reads with the best alignment of 2 or more mismatches and not parse reads that have better alignments. This allows for a stepwise examination of the taxonomic profiles generated by allowing gradually more mismatches. 
```
ngsLCA/ngsLCA -editdistmin 2 -editdistmax 2 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile.ed2
```
The similarity score function works exactly the same way as for the edit distance, with the exception that the similarity to the reference is now calculated as a percentage and hence takes the length of the read into account. Example for running the ngsLCA main program with similarity scores (95-100% similarity) [-simscore[low/high]] 0.95-1.0. 
```
ngsLCA/ngsLCA -simscorelow 0.95 -simscorehigh 1.0 -names ncbi_tax_dmp/names.dmp -nodes ncbi_tax_dmp/nodes.dmp -acc2tax ncbi_tax_dmp/nucl_gb.accession2taxid -bam sorted_bam_file.bam -outnames outfile.ss095to1
```

## The .lca output file format
The resulting file (.lca) from the ngsLCA main program is a regular text file in which each line contains a unique read that aligned to a reference with the given simlarity speficified when executing ngsLCA. 

The file is tab seperated and the first column contains the read metadata which is subdivided in a colon seperated style and the first is the read ID (colon seperated column 1-7), the query sequence (csc 8), the lenght of the sequence (csc 9), the number hits to a reference in the database(s) (scs 10). 

The second column contain the lowest taxonomic node the read has been assigned to, seperated by colon is the 'NCBI taxID':'taxon name':'taxonomic level' assigned to. Following columns contain the taxonomic path higher in the NCBI taxonomy for each assignment. 


# Running ngsLCA R package
## Example workflow
Running ngsLCA R package (and this example workflow) requires all lca files to be processed being placed in a folder as input, with each input lca file normally representing a sample. This example workflow will run through all ngsLCA R functions. The output files and results in this example will be generated into “Desktop/project_files/lca_files/run03”.

Define work directory and name of the run:

``` r
library(ngsLCA)

working_directory="Desktop/project_files/lca_files/" #change it to the path of your lca files folder.
run_name="run03"
```

Generate a complete taxa profile by combining all input lca files (this step is required for running all other functions):

``` r
ngsLCA_profile(path=working_directory,
               run=run_name,
               metadata="Desktop/project_files/sample_metadata.txt") #change it to the full path of your sample metadata.
```

Filter the combined taxa profile:

``` r
ngsLCA_filter(path=working_directory,
              run=run_name,
              remove.taxa="1115501,10114", #taxa ID indicate taxa that will be removed from the results.
              threshold.1=2,
              threshold.2=0,
              threshold.3=5)
```

De-contaminate the combined taxa profile:

``` r
ngsLCA_deContam(path=working_directory,
                run=run_name,
                control_path="Desktop/project_files/blank_control_lca_files/", #change it to the full path of your blank control lca files folder.
                threshold.1_control=2,
                threshold.2_control=0,
                threshold.3_control=5)
```

Group taxa into user-defined taxa groups:

``` r
ngsLCA_group(path=working_directory,
             run=run_name,
             group.name="Viruses,Archaea,Bacteria,Fungi,Viridiplantae,Metazoa",
             threshold.perGroup=0)
```

Classify taxa to taxonomic ranks:

``` r
ngsLCA_rank(path=working_directory,
            run=run_name,
            rank.name="species,genus,family")
```

Count reads and taxa numbers:

``` r
ngsLCA_count(path=working_directory,
             run=run_name)
```

Generate files that can be inputted into Megan for taxa profiling:

``` r
ngsLCA_meganFile(path=working_directory,
                 run=run_name)
```

Generate taxa abundance heatmaps:

``` r
ngsLCA_heatmap(path=working_directory,
               run=run_name,
               taxa.number=30)
```

Generate taxa abundance barplots:

``` r
ngsLCA_barplot(path=working_directory,
               run=run_name,
               taxa.number=20)
```

Perform rarefaction on taxa reads abundance:

``` r
ngsLCA_rarefy(path=working_directory,
              run=run_name)
```

Perform taxa profiles NMDS:

``` r
ngsLCA_NMDS(path=working_directory,
            run=run_name,
            dimension=2,
            trymax=10)
```

Generate taxa abundance stratplots:

``` r
ngsLCA_stratplot(path=working_directory,
                 run=run_name,
                 taxa.number=20)
```
