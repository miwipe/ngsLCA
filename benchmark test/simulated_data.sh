########################
########################
#simulate a dataset to test the parameters of ngsLCA
########################
########################


########################
#prepare genomes
########################
#download the assembly summary file from NCBI genome refseq manually

#make random number matix
shuf -i 1-246831 -n 20 > b_list.txt
shuf -i 1-444 -n 5 > f_list.txt
shuf -i 1-292 -n 3 > i_list.txt
shuf -i 1-148 -n 10 > p_list.txt
shuf -i 1-188 -n 5 > v_list.txt

#subset info from the refseq summary list
for file in *summary.txt
do
  bname=`echo $file | cut -d"_" -f1`
  awk 'NR == FNR{a[$0]; next};FNR in a' ${bname}_list.txt $file >> genome_list.txt
done

cut -f20 genome_list.txt > ftp_list.txt

#download genomes
cat ftp_list.txt | while read line
do
  ftp=`echo $line | sed "s/https/ftp/"`
  acc=`echo $ftp | cut -d\\/ -f10`
  echo ${ftp}/${acc}_genomic.fna.gz
  wget ${ftp}/${acc}_genomic.fna.gz
done

gzip -d *.gz



########################
#simulate random mutations with Mutation-Simulator
########################
for file in *.fna
do
  bname=`echo $file | cut -d. -f1`
  mutation-simulator $file args -sn 0.001 -titv 1
done



########################
#simulating
########################

#length distrubution
echo -e "30\t0.036257742
31\t0.03577441
32\t0.035312843
33\t0.0350033
34\t0.034423739
35\t0.033830851
36\t0.033269169
37\t0.032283047
38\t0.031280805
39\t0.030352632
40\t0.029164853
41\t0.028275628
42\t0.027571727
43\t0.026551788
44\t0.025688638
45\t0.025025628
46\t0.02407575
47\t0.023183003
48\t0.022431595
49\t0.021397358
50\t0.020590398
51\t0.020009471
52\t0.019326243
53\t0.018736239
54\t0.018259524
55\t0.017533403
56\t0.017013734
57\t0.016553989
58\t0.015857161
59\t0.015258202
60\t0.014792021
61\t0.01417582
62\t0.013702718
63\t0.013383309
64\t0.012986978
65\t0.012518551
66\t0.012217964
67\t0.011828796
68\t0.011431889
69\t0.011099245
70\t0.010643051
71\t0.010351206
72\t0.010086621
73\t0.009739619
74\t0.009460008
75\t0.009191112
76\t0.008915568
77\t0.008648251
78\t0.008439066
79\t0.008180644
80\t0.007914693" > len_dis.txt

#human contanimation
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gzip -d GCF_000001405.40_GRCh38.p14_genomic.fna.gz

#simulation
for file in *.fna
do
  bname=`basename $file | cut -d. -f1`
  mkdir $bname

  mkdir $bname/bact
  cp $path_to_gargammel/gargammel/bactDBexample/clovis/fasta/* $bname/bact/

  mkdir $bname/cont
  cp GCF_000001405.40_GRCh38.p14_genomic.fna $bname/cont/

  mkdir $bname/endo
  mv $file $bname/endo/

  gargammel -n 10000 -se --comp 0.08,0.02,0.9 -damagee 0.03,0.4,0.01,0.3 -damageb 0,0,0,0 -damagec 0,0,0,0 -qs 30 -f len_dis.txt -o $bname/simulation $bname

  cd $bname
  gzip -d simulation_d.fa.gz
  $path_to_BBMap/bbmap/reformat.sh in=simulation_d.fa out1=$bname.simulation.fq qfake=40
  mv $bname.simulation.fq ../
  mv endo/$file ../
  cd ../
  rm -r $bname
done



########################
#mapping simulated reads
########################

#mapping
mkdir mapping

cat *.simulation.fq > mapping/combined.simulation.fq

cd mapping
file=combined.simulation.fq
bname=$(basename "$file" | cut -d. -f1)
mkdir $bname
cd $bname

cat $NCBI_nt_and_refseq_database_list.txt | while read DB
do
bDB=$(basename $DB)
bowtie2 --threads 40 -k 1000 -x $DB -U ../$file --no-unal -S ${bname}.$bDB.sam 2> ${bname}.$bDB.log.txt
samtools view -@ 30 -b -o ${bname}.$bDB.sam.bam ${bname}.$bDB.sam
done

#sort bam files
samtools merge -n -@ 30 $bname.all_db_merged *.bam
mkdir TMP
samtools sort -n -T TMP -@ 30 -m 4G -O bam -o $bname.all_db_merged.sorted.bam $bname.all_db_merged
rm -r TMP


########################
#random subset reads to ngsLCA
########################

#seprate bam files
samtools view $bname.all_db_merged.sorted.bam -H > combined.bam_header.txt
samtools view $bname.all_db_merged.sorted.bam > combined.bam_alignment.txt

#subset simulated read headers
for file in *.simulation.fq
do
  bname=`basename $file | cut -d. -f1`
  grep "@" $file | cut -d"@" -f2 > $bname.readHeader
done


#simulate the dataset

run=run01
mkdir $run
for file in *.readHeader
do
  bname=`echo $file | cut -d. -f1`
  NO=`shuf -i 1-500 -n 1`
  shuf -i 1-10000 -n $NO > $bname.$run.list.txt
  awk 'NR == FNR{a[$0]; next};FNR in a' $bname.$run.list.txt $file > $bname.readHeader.txt
  rm $bname.$run.list.txt
  mv $bname.readHeader.txt $run/

  acc=`grep $bname $path_to_genome_acc2taxaID/genome_acc2taxaID.txt | cut -f1`
  taxaID=`grep $bname $path_to_genome_acc2taxaID/genome_acc2taxaID.txt | cut -f2`
  echo "$acc,$taxaID,$NO" >> $run/readsCompasation.csv
done

cd $run
cat *readHeader.txt > simulated_data.$run.txt

#subset the bam file
grep -f simulated_data.$run.txt ../combined.bam_alignment.txt > alignment.txt
cat ../combined.bam_header.txt alignment.txt | samtools view -b -o simulated_data.$run.bam


#ngsLCA
$path_to_ngsLCA/ngsLCA \
-names $path_to_NCBI_taxonomy/names.dmp \
-nodes $path_to_NCBI_taxonomy/nodes.dmp \
-acc2tax $path_to_genome_acc2taxaID/nucl_gb.accession2taxid \
-editdistmin 0 -editdistmax 0 -bam simulated_data.$run.bam -outnames simulated_data.$run.min0_max0

$path_to_ngsLCA/ngsLCA \
-names $path_to_NCBI_taxonomy/names.dmp \
-nodes $path_to_NCBI_taxonomy/nodes.dmp \
-acc2tax $path_to_genome_acc2taxaID/nucl_gb.accession2taxid \
-editdistmin 0 -editdistmax 1 -bam simulated_data.$run.bam -outnames simulated_data.$run.min0_max1

$path_to_ngsLCA/ngsLCA \
-names $path_to_NCBI_taxonomy/names.dmp \
-nodes $path_to_NCBI_taxonomy/nodes.dmp \
-acc2tax $path_to_genome_acc2taxaID/nucl_gb.accession2taxid \
-editdistmin 0 -editdistmax 2 -bam simulated_data.$run.bam -outnames simulated_data.$run.min0_max2
