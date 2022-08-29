#!/bin/bash
#test script for metadamage, 
#

PRG=../ngsLCA
BAM1=./data/f570b1db7c.dedup.filtered.bam


echo "--------------------"
echo "Using PRG: '${PRG}' and BAMFILE: '${BAM1}'"
echo "--------------------"

rm -f *.bin #remove old logfile and binary tempfile

RVAL=0

echo "Testing Existence of ${PRG}"
if [[ ! -f "${PRG}" ]]; then
    echo "Problem finding program: ${PRG}"
    RVAL=1
fi

echo "Testing Existence of ${BAM1}"
if [[ ! -f "${BAM1}" ]]; then
    echo "Problem finding file: ${PRG}"
    RVAL=$((2+RVAL))
fi

echo "Sorting bamfile"
BAM="$(dirname ${BAM1})/$(basename ${BAM1} .bam).rname.bam"
CMD="samtools sort -n ${BAM1} -o ${BAM}"
${CMD}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((4+RVAL))
fi

mkdir -p output

echo "Validating lca"
CMD="${PRG} -bam ${BAM} -names data/names.dmp.gz -nodes data/nodes.dmp.gz -acc2tax data/acc2taxid.map.gz -simscorelow 0.9 -simscorehigh 1.0 -editdistmin 0 -editdistmax 10 -minmapq 0 -fix_ncbi 0 -out output/test"
echo $CMD
${CMD} 
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((8+RVAL))
fi

echo "Validating checksum"
echo "========================"

sed 1d output/test.lca |md5sum -c files.md5
if [[ $? -ne 0 ]]; then
    echo "Problem with md5sum for lca file"
    RVAL=$((1024+${RVAL}))
fi

echo "=====RVAL:${RVAL}======="


if [[ ${RVAL} -ne 0 ]];then
    exit 1 #exit codes are capped at 255
fi
echo "All is ok, nothing failed"
exit 0
