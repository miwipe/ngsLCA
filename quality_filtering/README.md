# Raw read trimming and quality check 

This page decribes the process of trimming adaptors, quality filtering, removal of low complexity and deduplication of raw fastq data

We here use two tools (fastp and sga) for processing raw fastq data, however there exists a lot more out there and the reason for choosing these are that they are easy to use, fast, reliable and can perform the desired analysis or filtering. 
Both tools are available through bioconda, and we therefore recommend setting up conda and running the QC filtering in this environment. 

First we create the conda environment, and activate it
```
conda env create --file environment.yaml
conda activate ngsLCA
```

Next, we first use fastp (https://github.com/OpenGene/fastp) to trim the reads, remove poly X tails, low complexity reads and nucleotides with low qualities (OBS the deduplication function in version 0.23.1 isn't removing duplicates and we therefore added the steop below). 
```
for file in *fq
do
fastp -i $file -o $file.out.fq -V -D --dup_calc_accuracy 5  -g -x -q 30 -e 25 -l 30 -y -c -p -h $file.fastp.report.html -w 8
done 
```

Then sga (https://github.com/jts/sga) is used to remove duplicates reference free, defining duplicates strictly as reads that are identical or contained within another read. 
```
for file in *.out.fq
do
sga index --algorithm=ropebwt --threads=30 $file
sga filter --threads=30  --no-kmer-check $file -o $file.dedup.fq
done
```

The fastq files are now ready for downstream mapping and analysis. 

