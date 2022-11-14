# NJU-seq in *Arabidopsis thaliana*

## 1 Reference and index
### 1.1 *Arabidopsis thaliana* reference genome
+ Download reference
```bash
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O data/ath_dna.fa.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz -O data/ath.gff3.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz -O data/ath_transcript.fa.gz
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz -O data/ath_ncrna.fa.gz
```
+ Build index
```bash
# rRNA
cat NJU_seq/data/ath_rrna/* >data/ath_rrna.fa
bowtie2-build data/ath_rrna.fa index/ath_rrna

# mRNA
pigz -dc data/ath_transcript.fa.gz |
  perl NJU_seq/tool/fetch_fasta.pl \
  --stdin -s 'transcript_biotype:protein_coding' \
  >data/ath_protein_coding.fa
bowtie2-build data/ath_protein_coding.fa index/ath_protein_coding
```

### 1.2 Bacteria representative genome
+ Download
```bash
mkdir bacteria
cd bacteria

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ./representative.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/rsync.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
bash ASSEMBLY/representative.assembly.rsync.sh
bash ASSEMBLY/representative.assembly.collect.sh

# md5
cat ASSEMBLY/rsync.tsv |
    tsv-select -f 1 |
    parallel -j 4 --keep-order '
        echo "==> {}"
        cd ASSEMBLY/{}
        md5sum --check md5checksums.txt
    ' |
    grep -v ": OK"
```
+ Build index
```bash
# tRNA and rRNA and ncRNA(*rna_from_genomic.fna.gz)
JOB=$(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d)
for i in $JOB;do
  echo "=====>$i"
  gzip -ckd $i/*_rna_from_genomic.fna.gz >> bacteria_RNA.fa
done

# rRNA
for i in $JOB;do
  echo "=====>$i"
  gzip -ckd $i/*_rna_from_genomic.fna.gz |
    perl ../NJU_seq/tool/fetch_fasta.pl --stdin -s 'rRNA' >> bacteria_rRNA.fa
done

# quality_control(去除含有N的序列)
perl script/remove_sequence_with_N.pl bacteria_RNA.fa > tem&&
  mv tem bacteria_RNA.fa

perl script/remove_sequence_with_N.pl bacteria_rRNA.fa > tem&&
  mv tem bacteria_rRNA.fa

# build index
makie index
bowtie2-build ./bacteria_RNA.fa index/bacteria_RNA
bowtie2-build ./bacteria_rRNA.fa index/bacteria_rRNA
```

## 2 Data Selection and quality overview
+ Cut adapter and turn pair-end sequence data to single-end data
```bash
# Select data for analysing
for PREFIX in Ath_root_RF_NC Ath_root_RF_1 Ath_root_RF_2 Ath_root_RF_3;do
mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}" "raw_data/${PREFIX}"

# cut adapter
cutadapt -O 6 -m 10 -e 0.1 --discard-untrimmed -j 16 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-o data/${PREFIX}/R1.fq.gz \
-p data/${PREFIX}/R2.fq.gz \
raw_data/${PREFIX}/R1.fq.gz \
raw_data/${PREFIX}/R2.fq.gz

# Quality control for clean data
# Turn pair-end sequence data to single-end file
perl NJU_seq/quality_control/pe_consistency.pl \
  data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R2.fq.gz \
  temp/"${PREFIX}".fq.gz
done

#time perl NJU_seq/quality_control/fastq_qc.pl \
#  temp/Ath_root_RF_NC.fastq.gz \
#  temp/Ath_root_RF_1.fastq.gz \
#  temp/Ath_root_RF_2.fastq.gz \
#  temp/Ath_root_RF_3.fastq.gz \
#  output \
#  Ath_root_RF
```
+ remove bacteria genome
```bash
for TISSUE in flower leaf root stem;do
  for PREFIX in Ath_${TISSUE}_NC Ath_${TISSUE}_1 Ath_${TISSUE}_2 Ath_${TISSUE}_3;do
    mkdir -p output/${PREFIX}
    THREAD=24

    # align
    time bowtie2 -p "${THREAD}" -a -t \
      --end-to-end -D 20 -R 3 \
      -N 0 -L 10 -i S,1,0.50 --np 0 \
      --xeq -x index/bacteria_rRNA \
      -1 ../data/"${PREFIX}"/R1.fq.gz -2 ../data/"${PREFIX}"/R2.fq.gz \
      -S output2/"${PREFIX}"/bacteria_align.sam \
      2>&1 |
      tee output2/"${PREFIX}"/bacteria.bowtie2.log

    time pigz -p "${THREAD}" output/"${PREFIX}"/bacteria_align.sam
  done
done

# remove
for TISSUE in flower leaf root stem;do
    for PREFIX in Ath_${TISSUE}_NC Ath_${TISSUE}_1 Ath_${TISSUE}_2 Ath_${TISSUE}_3;do
      THREAD=24

      pigz -dcf output/"${PREFIX}"/bacteria_align.sam.gz  |
        grep -v "@" |
        parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
          awk '\''$6!="*"&&$7=="="{print $1 "\t" $6}
          '\'' |
          perl script/align_filter.pl 
        ' | uniq > output/${PREFIX}/filter_name.tsv


  done
done

find job_list -maxdepth 1 -mindepth 1 -type f |
  parallel --linebuffer -k -j 4 "
    echo >&2 '====> {1}'
    perl script/remove_select_reads.pl -r {1} -f output/${PREFIX}/filter_name.tsv
  " > output/${PREFIX}/test.fastq
```


## 3 Alignment and Filter
```bash
THREAD=16
PREFIX='Ath_root_RF_NC'

for PREFIX in Ath_root_RF_NC Ath_root_RF_1 Ath_root_RF_2 Ath_root_RF_3;do
time bowtie2 -p "${THREAD}" -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz \
  -S output/"${PREFIX}"/rrna.raw.sam \
  2>&1 |
  tee output/"${PREFIX}"/rrna.bowtie2.log

perl NJU_seq/tool/stat_alignment.pl \
  output/"${PREFIX}"/rrna.bowtie2.log |
  Rscript NJU_seq/tool/draw_table.R \
  output/"${PREFIX}"/rrna.bowtie2.pdf

time pigz -p "${THREAD}" output/"${PREFIX}"/rrna.raw.sam
done  
```

+ Filter
```bash
THREAD=16
PREFIX='Ath_root_RF_NC'

time pigz -dcf output/"${PREFIX}"/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
  >temp/"${PREFIX}"/rrna.out.tmp
# matchquality_judge.pl ：
# LENGTH < 20 and all matched、LENGTH is [20,30) and at most 1 mismatch、LENGTH is [30,40) and at most 2 mismatches、LENGTH > 40 and at most 3 mismatches

# multimatch_judge.pl ：
# 一个reads如果只有一种匹配,则直接输出,如果有多种匹配,有完全匹配的就输出完全匹配,没有完全匹配的就输出匹配长度最长的那种
```

