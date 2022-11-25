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

# quality_control(去除含有N的序列)
perl script/remove_sequence_with_N.pl bacteria_RNA.fa > tem&&
  mv tem bacteria_RNA.fa

# build index
makie index
bowtie2-build --threads 8 ./bacteria_RNA.fa index/bacteria_RNA
rm bacteria_RNA.fa
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
```

+ raw sequence file statistics
```bash
JOB=$(find data -maxdepth 1 -mindepth 1 -type d)

for J in $JOB;do
  echo "===> $J"

  pigz -dcf $J/R1.fq.gz | perl -e'
  while (<>) {
    chomp( my $seq_name = $_ );
    chomp( my $seq = <> );
    my $seq_length = length($seq);
    chomp( my $info = <> );
    chomp( my $quality = <> );
    print "$seq_length\n";
  }
' | sort -n | uniq -c | perl -ne'
      s/^\s+//;
      print "$_";
    ' | tr " " "\t" | tsv-select --fields 2,1 > $J/length_distribution.tsv
  
  (echo -e "reads_length\tTotal_number" && cat $J/length_distribution.tsv) > tem&&
    mv tem $J/length_distribution.tsv
done
```

### 2.1 1331种物种代表性菌株一一匹配
```bash
bash script/bacteria_align.sh
bash script/statistics.sh
bash script/merge_plot.sh

# sort and merge
cd bacteria
cat output/${PREFIX}/align_statistics.tsv | keep-header -- sort -nr -k4,4 > tem&&
  mv tem output/${PREFIX}/align_statistics.tsv

tsv-join --filter-file <(cut -d "," -f 1,2 Bacteria.assembly.collect.csv | tr "," "\t") -H --key-fields name --append-fields Organism_name output/${PREFIX}/align_statistics.tsv | 
tsv-select -H -f Organism_name --rest last > tem&&
  mv tem output/${PREFIX}/align_statistics.tsv

cat output/${PREFIX}/$J/reads_info.tsv >> output/${PREFIX}/total_remove_reads_info.tsv
cat output/${PREFIX}/total_remove_reads_info.tsv | sort | uniq > tem&&
  mv tem output/${PREFIX}/total_remove_reads_info.tsv

cut -f 3 output/${PREFIX}/total_remove_reads_info.tsv | 
  sort -n | uniq -c |
  perl -ne'
    s/^\s+//;
    print "$_";
  ' |
  tr " " "\t" | tsv-select --fields 2,1 > output/${PREFIX}/length_distribution.tsv

(echo -e "reads_length\tRemove_number" && cat output/${PREFIX}/length_distribution.tsv) > tem&&
  mv tem output/${PREFIX}/length_distribution.tsv 

# plot
cd NJU_seq_analysis_ath
tsv-join -H --filter-file bacteria/output/${PREFIX}/length_distribution.tsv --key-fields reads_length --append-fields Remove_number data/${PREFIX}/length_distribution.tsv

perl bacteria/script/tsv_join_plus.pl data/${PREFIX}/length_distribution.tsv bacteria/output/${PREFIX}/length_distribution.tsv > tem&&
  mv tem data/${PREFIX}/length_distribution.tsv

echo -e "reads_length\tnumber\tgroup" > data/${PREFIX}/plot.tsv
sed '1d' data/${PREFIX}/length_distribution.tsv | perl -ne'
  chomp;
  if (/^(\S+)\t(\S+)\t(\S+)/) {
    my $reads_length = $1;
    my $total_number = $2;
    my $remove_number = $3;
    my $keep_number = ( $total_number - $remove_number );
    print "$reads_length\t$keep_number\tKeep\n";
    print "$reads_length\t$remove_number\tRemove\n";
  }
' >> data/${PREFIX}/plot.tsv
```
```R
library("ggplot2")
DATA <- read.table("plot.tsv", header = TRUE, sep = "\t")

ggplot(DATA, aes(reads_length, number, fill=group)) +
  geom_bar(stat = "identity")+
  xlim(9,60)+
  theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold"))
```

+ remove
```bash
pigz -dcf R1.fq.gz | grep "@" | grep -v -w -f discard.tsv > keep.tsv
pigz -dcf R1.fq.gz | 
  grep -A3 -w -f keep.tsv | 
  perl -ne'
    if ( ! (/^--$/) ) {
      print "$_";
    }
  ' > R1_filter.fq

pigz -p "${THREAD}" R1_filter.fq
```

### 2.2 与rRNA恒定区匹配
**细菌的16SrDNA中有多个区段保守性，根据这些保守区可以设计出细菌通用物，可以扩增出所有细菌的16SrDNA片段，并且这些引物仅对细菌是特异性的，也就是说这些引物不会与非细菌的DNA互补**
```bash
mkdir -p silva



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
