########################
########################
#compare ngsLCA to the other 2 available programs
########################
########################

mkdir bam
#download the bam files from European Nucleotide Archive under accession number PRJEB14494
cd ..


########################
#test sam2rma
########################

mkdir sam2rma
cd sam2rma
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip


#mornitoring memory usage
while true
do
  ps wwaux | grep "java" | head -n1 >> SPL.sam2rma_memmory_log.txt
  sleep 2
done

#run sam2rma
for file in ../bam/*.bam
do
  bname=`basename $file | cut -d. -f1`
  echo processing $file
  samtools view $file -h > $bname.sam
  date
  $path_to_megan/tools/sam2rma -i $bname.sam -o $bname.sam2rma -a2t nucl_acc2tax-Jul2019.abin \
  -ms 0 -supp 0 -sup 1 -t 1
  date
  rm $bname.sam
done &>sam2rma.log.txt

cd ..


########################
#test sam2lca
########################

#index the bam files
cd bam
for file in *.bam
do
  bname=`basename $file | cut -d. -f1`
  samtools sort -@ 20 -o $bname.q_sorted.bam $file
  samtools index -@ 20 -c $bname.q_sorted.bam
done

#run sam2lca
mkdir sam2lca
cd sam2lca

#mornitoring memory usage
while true
do
  ps wwaux | grep "sam2lca" | head -n1 >> SPL.sam2lca_memmory_log.txt
  sleep 2
done

for file in ../bam/*.q_sorted.bam
do
  echo start $file
  date
  sam2lca analyze -a nucl -i 1 -p 1 $file
  date
  echo end $file
done &>SPL.sam2lca_time_log.txt


########################
#test ngsLCA
########################

mkdir ngsLCA
cd ngsLCA

#mornitoring memory usage
while true
do
  ps wwaux | grep "ngsLCA" | head -n1 >> SPL.ngsLCA_memmory_log.txt
  sleep 2
done

for file in ../bam/*.bam
do
  processing $file
  date
  $path_to_ngsLCA/ngsLCA \
  -names $path_to_NCBI_taxonomy/names.dmp \
  -nodes $path_to_NCBI_taxonomy/db_files/nodes.dmp \
  -acc2tax $path_to_NCBI_acc2taxID/nucl_gb.accession2taxid \
  -editdistmin 0 -editdistmax 0 -bam $file -outnames $file
  date
done &>SPL.ngsLCA_time_log.txt
