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
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 2 |
  parallel -j 24 --keep-order "
    echo '===>index_{1}'
    makdir -p index/{1}
    bowtie2-build ASSEMBLY/{1}/*_rna_from_genomic.fna.gz index/{1}/{1} 
  "
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

JOB=$(ls output/${PREFIX})
for J in $JOB;do
  echo "===> $J"
  cat output/${PREFIX}/$J/reads_info.tsv >> output/${PREFIX}/total_remove_reads_info.tsv
done

tsv-uniq -f 1 output/${PREFIX}/total_remove_reads_info.tsv > tem&&
  mv tem output/${PREFIX}/total_remove_reads_info.tsv

cut -f 3 output/${PREFIX}/total_remove_reads_info.tsv | 
  sed '1d' |
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
  xlim(9,70)+
  labs(title = "Ath_flower_NC")+
  theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold"), plot.title = element_text(hjust=0.5))
```

+ remove
```bash
cd NJU_seq_analysis_ath
PREFIX=Ath_flower_1
cut -f 1 bacteria/output/${PREFIX}/total_remove_reads_info.tsv | sed '1d' > data/${PREFIX}/discard.lst
pigz -dcf data/${PREFIX}/R1.fq.gz | grep "@" |
  cut -d " " -f 1 | sed 's/^@//g' |
  grep -v -w -f data/${PREFIX}/discard.lst > data/${PREFIX}/keep_ID.lst

brew install seqkit
for J in R1 R2;do
  echo "====> $J"
  seqkit grep -f data/${PREFIX}/keep_ID.lst data/${PREFIX}/${J}.fq.gz > data/${PREFIX}/${J}_filter.fq
  pigz -p 4 data/${PREFIX}/${J}_filter.fq
done
```

### 2.2 与rRNA恒定区匹配
**细菌的16SrDNA中有多个区段保守性，根据这些保守区可以设计出细菌通用物，可以扩增出所有细菌的16SrDNA片段，并且这些引物仅对细菌是特异性的，也就是说这些引物不会与非细菌的DNA互补**
```bash
cd bacteria
mkdir -p rRNA_conserve

# 提取所有菌株的16SrRNA与23SrRNA
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 2 |
  parallel -j 2 --keep-order '
    echo "===> extract_{1}"

    pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | faops some stdin <(pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | grep ">" | grep "16S" | cut -d ">" -f 2 | head -n 1) stdout >> rRNA_conserve/16S.fas

    pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | faops some stdin <(pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | grep ">" | grep "23S" | cut -d ">" -f 2 | head -n 1) stdout >> rRNA_conserve/23S.fas
  '

perl script/remove_sequence_with_N.pl rRNA_conserve/16S.fas > tem&&
  mv tem rRNA_conserve/16S.fas
perl script/remove_sequence_with_N.pl rRNA_conserve/23S.fas > tem&&
  mv tem rRNA_conserve/23S.fas

# 多序列比对
brew install mafft

mafft --reorder rRNA_conserve/16S.fas > rRNA_conserve/16S_align.fas
mafft --reorder rRNA_conserve/23S.fas > rRNA_conserve/23S_align.fas

# 保守区查找

```

### 2.3 不同环境细菌的mRNA匹配情况
```bash
mlr --itsv --omd cat
```

| **植物相关** |  |  |  |
| --- | --- | --- | --- |
| Lactococcus lactis | 乳酸乳球菌 | 最主要存在植物表面存在，被摄入人体内之后才在肠道定植 | 
| Xanthomonas campestris | 野油菜黄单胞菌 | 可引起多种植物疾病的细菌，十字花科黑腐病等 |
| Streptococcus thermophilus | 嗜热链球菌 | 菌株能大量的从植物中分离，并分布到环境中，其降解特性可用于乳制品 |
| Pseudomonas syringae | 丁香假单胞菌 | 植物病原体 |
| Pseudomonas protegens | 恶臭假单胞菌 | 植物保护细菌，能够产生对抗植物病原体的物质 |
| Bacillus velezensis | 贝雷兹芽孢杆菌 | 保护植物不受病原菌侵袭，植物内生菌 |
| Serratia marcescens | 粘质沙雷氏菌 | 从多种植物中都能分离得到 |
| **水体土壤等环境** |  |  |
| Clostridium acetobutylicum | 丙酮丁醇梭菌 | 最常生活在土壤中，不过也在多种环境中被鉴定到 |
| Rhodococcus erythropolis | 红球菌 | 通常在土壤和水体中发现，以及真核细胞，野生环境中居多，利用多种有机物 |
| Arcobacter butzleri | 布氏弧杆菌 | 该细菌在环境中广泛存在，主要是各种水体分离得到该菌株，会引起腹泻 |
| Clostridium beijerinckii | 拜氏梭菌 | 在自然界中普遍存在，并且最初分离自土壤样品 |
| Azotobacter vinelandii | 维氏固氮杆菌 | 革兰氏阴性的固氮细菌，严格需氧且广泛自由生存的土壤微生物 |
| Acinetobacter baumannii | 鲍曼不动杆菌 | 其大量存在于环境中，在水体，土壤甚至植物表面都有发现 |
| **人体口腔、皮肤等** |  |  |
| Escherichia coli | 大肠杆菌 | 最主要是肠道定植，也广泛存在与环境中 |
| Streptococcus mutans | 变形链球菌 | 引起龋齿，主要生活在人类口腔 |
| Staphylococcus epidermidis | 表皮葡萄球菌 | 每个人表皮都有的一种机会致病菌 |
| Staphylococcus aureus | 金黄色葡萄球菌 | 自然栖息地在人类和动物中，是皮肤自然菌群的一部分，当然也在环境中能够发现其存在 |
| Bifidobacterium breve  | 短双歧杆菌 | 主要生活在人类的及与人有关的生境中 |
| Streptococcus salivarius | 唾液链球菌 | 存在于人体上呼吸道和口腔 |
| **其他** |  |  |
| Riemerella anatipestifer | 鸭瘟立默氏菌 | 主要存在于鸭中，并且能够感染鸭 |
| Methanosarcina barkeri | 巴氏甲烷八叠球菌 | 在淡水湖Fusaro lake中的泥浆样本中发现了该菌株，能够产甲烷，是该属最常见的菌 |
| Methanosarcina siciliae | 西西里甲烷八叠球菌 | 从海底泥中分离的厌氧菌 |

```bash
cd bacteria_mRNA

# index
mkdir index

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 2 |
  parallel -j 2 --keep-order '
    echo "===> index_{1}"
    mkdir -p index/{1}

    bowtie2-build --threads 2 ASSEMBLY/{1}/*_cds_from_genomic.fna.gz index/{1}/{1}
  '

# align
mkdir output

for PREFIX in Ath_root_NC;do
  find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    cut -d "/" -f 2 |
    parallel -j 2 --line-buffer "
      echo '====> align_{1}'
      mkdir -p output/${PREFIX}/{1}

      bowtie2 -a -t \
      --end-to-end -D 20 -R 3 \
      -N 0 -L 10 -i S,1,0.50 --np 0 \
      --xeq -x index/{1}/{1} \
      -1 ../data/${PREFIX}/R1.fq.gz -2 ../data/${PREFIX}/R2.fq.gz \
      -S output/${PREFIX}/{1}/{1}_align.sam \
      2>&1 |
      tee output/${PREFIX}/{1}/{1}.bowtie2.log 

      pigz output/${PREFIX}/{1}/{1}_align.sam
    "
done

# statistics_mRNA
bash script/statistics.sh
bash script/merge_plot.sh

cat output/${PREFIX}/align_statistics.tsv | keep-header -- sort -nr -k4,4 > ${PREFIX}_mRNA.tsv

tsv-join --filter-file <(cut -d "," -f 1,2 ../bacteria/Bacteria.assembly.collect.csv | tr "," "\t") -H --key-fields name --append-fields Organism_name ${PREFIX}_mRNA.tsv | 
tsv-select -H -f Organism_name --rest last > tem&&
  mv tem ${PREFIX}_mRNA.tsv

cut -f 3 output/${PREFIX}/total_remove_reads_info.tsv | 
  sed '1d' |
  sort -n | uniq -c |
  perl -ne'
    s/^\s+//;
    print "$_";
  ' |
  tr " " "\t" | tsv-select --fields 2,1 > output/${PREFIX}/length_distribution.tsv

(echo -e "reads_length\tRemove_number" && cat output/${PREFIX}/length_distribution.tsv) > tem&&
  mv tem output/${PREFIX}/length_distribution.tsv 

perl ../bacteria/script/tsv_join_plus.pl data/${PREFIX}/length_distribution.tsv output/${PREFIX}/length_distribution.tsv > tem&&
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

# statistics_rRNA_tRNA
cd bacteria_mRNA

( head -n 1 ../bacteria/output/${PREFIX}/align_statistics.tsv && cat ../bacteria/output/${PREFIX}/align_statistics.tsv | grep -w -f <(ls ASSEMBLY) ) > tem&&
  mv tem output/${PREFIX}_rRNA_tRNA.tsv
```

## 3 Alignment and Filter
```bash
THREAD=24

for PREFIX in Ath_flower_NC Ath_flower_1 Ath_flower_2 Ath_flower_3;do
  echo "====> ${PREFIX}"
  mkdir -p output_filter/${PREFIX}

  bowtie2 -p ${THREAD} -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/${PREFIX}/R1_filter.fq.gz -2 data/${PREFIX}/R2_filter.fq.gz \
  -S output_filter/${PREFIX}/rrna.raw.sam \
  2>&1 |
  tee output_filter/${PREFIX}/rrna.bowtie2.log

  pigz -p ${THREAD} output_filter/${PREFIX}/rrna.raw.sam
done  
```

+ Filter
```bash
THREAD=4
PREFIX=Ath_root_NC

time pigz -dcf output_filter/${PREFIX}/rrna.raw.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j ${THREAD} '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
  >temp/${PREFIX}/filter_rrna.out.tmp
# matchquality_judge.pl ：
# LENGTH < 20 and all matched、LENGTH is [20,30) and at most 1 mismatch、LENGTH is [30,40) and at most 2 mismatches、LENGTH > 40 and at most 3 mismatches

# multimatch_judge.pl ：
# 一个reads如果只有一种匹配,则直接输出,如果有多种匹配,有完全匹配的就输出完全匹配,没有完全匹配的就输出匹配长度最长的那种

time parallel -j 3 "
  perl NJU_seq/rrna_analysis/readend_count.pl \\
    NJU_seq/data/ath_rrna/{}.fa temp/${PREFIX}/filter_rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s
```

```bash
perl ../bacteria/script/tsv_join_plus.pl data/${PREFIX}/length_distribution.tsv output/${PREFIX}/length_distribution.tsv > tem&&
  mv tem data/${PREFIX}/length_distribution.tsv
```