#! usr/bin/bash
# bash align_rRNA_mRNA.sh raw_output/stage_align/remove_align ${PREFIX} rRNA_conserve/rRNA_tRNA/mRNA

DIR=$1
PREFIX=$2
STAGE=$3

if [ $1 == stage_align ]; then
    DIR2=filter
elif [ $1 == remove_align ]; then
    DIR2=remove
fi

if [ $1 == raw_output ]; then
    pigz -dcf ${1}/${PREFIX}/rrna.raw.sam.gz | grep -v "@" | tsv-filter --str-ne 6:* | cut -f 1 | uniq > ${1}/${PREFIX}/rrna.lst
    pigz -dcf ${1}/${PREFIX}/mrna.raw.sam.gz | 
        parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
            grep -v "@" | tsv-filter --str-ne 6:* | cut -f 1 | uniq
        ' > ${1}/${PREFIX}/mrna.lst
    cat ${1}/${PREFIX}/mrna.lst | sort | uniq > ${1}/${PREFIX}/tem&&
        mv ${1}/${PREFIX}/tem ${1}/${PREFIX}/mrna.lst
    
    echo "Calculate"
    COMMOM_READS=$(perl remove_bacteria/script/common.pl --index ${1}/${PREFIX}/rrna.lst --search ${1}/${PREFIX}/mrna.lst | wc -l)
    TOTAL_READS=$(pigz -dcf data/${PREFIX}/R1.fq.gz | grep "@" | wc -l)

    RATE=$(echo -e "scale=2;(${COMMOM_READS} * 100 / ${TOTAL_READS})" | bc)
    echo -e "Common reads ratio:${RATE}%" > ${1}/${PREFIX}/common.tsv
else
    pigz -dcf ${1}/${PREFIX}/${3}/rrna_align.sam.gz | grep -v "@" | tsv-filter --str-ne 6:* | cut -f 1 | uniq > ${1}/${PREFIX}/${3}/rrna.lst
    pigz -dcf ${1}/${PREFIX}/${3}/mrna_align.sam.gz | 
        parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
            grep -v "@" | tsv-filter --str-ne 6:* | cut -f 1 | uniq
        ' > ${1}/${PREFIX}/${3}/mrna.lst
    cat ${1}/${PREFIX}/${3}/mrna.lst | sort | uniq > ${1}/${PREFIX}/${3}/tem&&
        mv ${1}/${PREFIX}/${3}/tem ${1}/${PREFIX}/${3}/mrna.lst
    
    echo "Calculate"
    COMMOM_READS=$(perl remove_bacteria/script/common.pl --index ${1}/${PREFIX}/${3}/rrna.lst --search ${1}/${PREFIX}/${3}/mrna.lst | wc -l)
    TOTAL_READS=$(pigz -dcf remove_bacteria/${3}/${DIR2}/${PREFIX}/R1*fq.gz | grep "@" | wc -l)

    RATE=$(echo -e "scale=2;(${COMMOM_READS} * 100 / ${TOTAL_READS})" | bc)
    echo -e "Common reads ratio:${RATE}%" > ${1}/${PREFIX}/${3}/common.tsv
fi
