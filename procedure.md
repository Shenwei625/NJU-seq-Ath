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
+ Info
```bash
# Bacteria
echo -e "name\tBacteria_CG_content_rRNA_tRNA\tBacteria_CG_content_mRNA" > ASSEMBLY/statistics.tsv

JOB=$(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | cut -d "/" -f 2 )

for J in $JOB;do
echo "===> $J"
    export i=$J

    pigz -dcf ASSEMBLY/$J/*_rna_from_genomic.fna.gz | faops count stdin | tail -n 1 |
      perl -a -F"\t" -e '
        print "$ENV{i}\t";
        my $CG_content = ((( @F[3] + @F[4] ) / @F[1] ) * 100 );
        printf "%.2f%%\t", $CG_content;
      ' >> ASSEMBLY/statistics.tsv

    pigz -dcf ASSEMBLY/$J/*_cds_from_genomic.fna.gz | faops count stdin | tail -n 1 |
      perl -a -F"\t" -e '
        my $CG_content = ((( @F[3] + @F[4] ) / @F[1] ) * 100 );
        printf "%.2f%%\n", $CG_content;
      ' >> ASSEMBLY/statistics.tsv  
done

# Ath
cd Ath

#rrna
cat data/ath_rrna.fa | faops count stdin | tail -n 1 |
  perl -a -F"\t" -e '
    my $CG_content = ((( @F[3] + @F[4] ) / @F[1] ) * 100 );
    printf "%.2f%%\n", $CG_content;
  '
# 53.45%

#mrna
cat data/ath_protein_coding.fa | faops count stdin | tail -n 1 |
  perl -a -F"\t" -e '
    my $CG_content = ((( @F[3] + @F[4] ) / @F[1] ) * 100 );
    printf "%.2f%%\n", $CG_content;
  '
# 41.54%
```

## 2 Remove rRNA conserve region
```bash
cd remove_bacteria
mkdir -p remove_bacteria/rRNA_conserve
# extract 16S and 23S 
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 2 |
  parallel -j 2 --keep-order '
    echo "===> extract_{1}"

    pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | faops some stdin <(pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | grep ">" | grep "16S" | cut -d ">" -f 2 | head -n 1) stdout >> rRNA_conserve/16S.fas

    pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | faops some stdin <(pigz -dcf ASSEMBLY/{1}/*_rna_from_genomic.fna.gz | grep ">" | grep "23S" | cut -d ">" -f 2 | head -n 1) stdout >> rRNA_conserve/23S.fas
  '

# Remove sequence with N
perl script/remove_sequence_with_N.pl rRNA_conserve/16S.fas > tem&&
  mv tem rRNA_conserve/16S.fas
perl script/remove_sequence_with_N.pl rRNA_conserve/23S.fas > tem&&
  mv tem rRNA_conserve/23S.fas

# search conserve region(MAST,dowload MAST text output,16S/23S_region.tsv)
# 16S
perl script/16S_motif.pl 16S_region.tsv | sed 's/,$//g' > rRNA_conserve/16S_filter.tsv
cat rRNA_conserve/16S_filter.tsv | wc -l | xargs seq |
  parallel -j 4 -k --linebuffer '
    faops region rRNA_conserve/16S.fa <(sed -n "{1}p" rRNA_conserve/16S_filter.tsv) stdout >> rRNA_conserve/16S_conserve.fa
  '
seqkit rmdup -s rRNA_conserve/16S_conserve.fa -o rRNA_conserve/16S_conserve_rmdup.fa
# 23S
perl script/23S_motif.pl rRNA_conserve/23S_region.tsv | sed 's/,$//g' > rRNA_conserve/23S_filter.tsv
cat rRNA_conserve/23S_filter.tsv | wc -l | xargs seq |
  parallel -j 4 -k --linebuffer '
    faops region rRNA_conserve/23S.fa <(sed -n "{1}p" rRNA_conserve/23S_filter.tsv) stdout >> rRNA_conserve/23S_conserve.fa
  '
seqkit rmdup -s rRNA_conserve/23S_conserve.fa -o rRNA_conserve/23S_conserve_rmdup.fa

# index
mkdir -p index/16S index/23S 
bowtie2-build --threads 2 rRNA_conserve/16S_conserve_rmdup.fa index/16S/16S
bowtie2-build --threads 2 rRNA_conserve/23S_conserve_rmdup.fa index/23S/23S

# align
mkdir -p rRNA_conserve/output
PREFIX=Ath_flower_NC
for J in 16S 23S;do
    mkdir -p rRNA_conserve/output/${PREFIX}/${J}

    bowtie2 -p 20 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/${J}/${J} \
    -1 ../data/${PREFIX}/R1.fq.gz -2 ../data/${PREFIX}/R2.fq.gz \
    -S rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam \
    2>&1 |
    tee rRNA_conserve/output/${PREFIX}/${J}/${J}.bowtie2.log 

    pigz -p 20 rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam
done

# statistics
for J in 16S 23S;do
  echo "====> $J"

  pigz -dcf rRNA_conserve/output/${PREFIX}/${J}/${J}_align.sam.gz |
    grep -v "@" |
    tsv-filter --regex '6:^[0-9]+=$' |
    tsv-filter --str-eq 7:= |
    cut -f 1,10 | tsv-uniq -f 1 > rRNA_conserve/output/${PREFIX}/${J}/${J}_match.tsv

  perl script/pre_remove_statistics.pl ../data/${PREFIX}/R1.fq.gz rRNA_conserve/output/${PREFIX}/${J}/${J}_match.tsv \
  rRNA_conserve/output/${PREFIX}/${J}/reads_info.tsv \
  rRNA_conserve/output/${PREFIX}/${J}/remove_info.tsv
done

# merge
bash script/merge.sh rRNA_conserve $PREFIX

# remove
mkdir -p rRNA_conserve/filter/${PREFIX}

pigz -dcf ../data/${PREFIX}/R1.fq.gz | grep "@" |
    cut -d " " -f 1 | sed 's/^@//g' |
    grep -v -w -f <(sed '1d' rRNA_conserve/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) > rRNA_conserve/filter/${PREFIX}/keep.lst

for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f rRNA_conserve/filter/${PREFIX}/keep.lst ../data/${PREFIX}/${J}.fq.gz > rRNA_conserve/filter/${PREFIX}/${J}_conserve.fq

  pigz -p 4 rRNA_conserve/filter/${PREFIX}/${J}_conserve.fq
done
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```

## 3 Remove rRNA and tRNA of 1331 bacteria
```bash
cd remove_bacteria
mkdir -p rRNA_tRNA
mkdir -p index/rRNA_tRNA

# index tRNA and rRNA and ncRNA(*_rna_from_genomic.fna.gz)
find ../../NJU_seq_analysis_ath/bacteria/ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 6 |
  parallel -j 24 --keep-order "
    echo '===>index_{1}'
    mkdir -p index/rRNA_tRNA/{1}
    bowtie2-build ../../NJU_seq_analysis_ath/bacteria/ASSEMBLY/{1}/*_rna_from_genomic.fna.gz index/rRNA_tRNA/{1}/{1} 
  "

# align
PREFIX=Ath_flower_NC
bash script/align_bacteria.sh ${PREFIX}

# statistics
bash script/statistics.sh rRNA_tRNA $PREFIX

# merge
bash script/merge.sh rRNA_tRNA $PREFIX

# remove
mkdir -p rRNA_tRNA/filter/${PREFIX}

pigz -dcf rRNA_conserve/filter/${PREFIX}/R1_conserve.fq.gz | grep "@" |
    cut -d " " -f 1 | sed 's/^@//g' |
    grep -v -w -f <(sed '1d' rRNA_tRNA/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) > rRNA_tRNA/filter/${PREFIX}/keep.lst

for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f rRNA_tRNA/filter/${PREFIX}/keep.lst rRNA_conserve/filter/${PREFIX}/${J}_conserve.fq.gz > rRNA_tRNA/filter/${PREFIX}/${J}_conserve_rRNA.fq

  pigz -p 4 rRNA_tRNA/filter/${PREFIX}/${J}_conserve_rRNA.fq
done
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```

## 4 Remove mRNA of 1331 bacteria
```bash
cd remove_bacteria
mkdir -p index/mRNA

# index mRNA(*_cds_from_genomic.fna.gz)
find ../../NJU_seq_analysis_ath/bacteria/ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 6 |
  parallel -j 2 --keep-order "
    echo '===>index_{1}'
    mkdir -p index/mRNA/{1}

    bowtie2-build --threads 2 ../../NJU_seq_analysis_ath/bacteria/ASSEMBLY/{1}/*_cds_from_genomic.fna.gz index/mRNA/{1}/{1} 
  "

# split_align
PREFIX=Ath_flower_NC
mkdir -p mRNA/output/${PREFIX}
mkdir -p job

find index/mRNA -maxdepth 1 -mindepth 1 -type d |
  cut -d "/" -f 3 |
  split -l 100 -a 2 -d - job/ 

JOB=$(find job -maxdepth 1 -type f -name "[0-9]?" | sort -n)
for J in $JOB;do
  bsub -q largemem -n 24 -J "$J" "bash script/align_mRNA.sh ${PREFIX} $J"
done

# statistics 
bash script/statistics.sh mRNA ${PREFIX}

# merge
bash script/merge.sh mRNA ${PREFIX}

# remove
mkdir -p mRNA/filter/${PREFIX}

pigz -dcf rRNA_tRNA/filter/${PREFIX}/R1_conserve_rRNA.fq.gz | grep "@" |
    cut -d " " -f 1 | sed 's/^@//g' |
    grep -v -w -f <(sed '1d' mRNA/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) > mRNA/filter/${PREFIX}/keep.lst

for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f mRNA/filter/${PREFIX}/keep.lst rRNA_tRNA/filter/${PREFIX}/${J}_conserve_rRNA.fq.gz > mRNA/filter/${PREFIX}/${J}_conserve_rRNA_mRNA.fq

  pigz -p 4 mRNA/filter/${PREFIX}/${J}_conserve_rRNA_mRNA.fq
done
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```

## 5 Align
+ Different remove stages
```bash
cd Ath
PREFIX=Ath_flower_NC
mkdir -p stage_align

bash remove_bacteria/script/align.sh stage rRNA_conserve ${PREFIX} rrna
bash remove_bacteria/script/align.sh stage rRNA_tRNA ${PREFIX}  rrna
bash remove_bacteria/script/align.sh stage mRNA ${PREFIX} rrna

bash remove_bacteria/script/align.sh stage rRNA_conserve ${PREFIX} protein_coding
bash remove_bacteria/script/align.sh stage rRNA_tRNA ${PREFIX}  protein_coding
bash remove_bacteria/script/align.sh stage mRNA ${PREFIX} protein_coding
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```

+ Remove reads
```bash
# extract
cd remove_bacteria
PREFIX=Ath_flower_NC

mkdir -p rRNA_conserve/remove/${PREFIX}
for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f <(sed '1d' rRNA_conserve/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1 ) ../data/${PREFIX}/${J}.fq.gz > rRNA_conserve/remove/${PREFIX}/${J}_conserve.re.fq

  pigz -p 4 rRNA_conserve/remove/${PREFIX}/${J}_conserve.re.fq
done

mkdir -p rRNA_tRNA/remove/${PREFIX}
for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f <(sed '1d' rRNA_tRNA/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1) rRNA_conserve/filter/${PREFIX}/${J}_conserve.fq.gz > rRNA_tRNA/remove/${PREFIX}/${J}_conserve_rRNA.re.fq

  pigz -p 4 rRNA_tRNA/remove/${PREFIX}/${J}_conserve_rRNA.re.fq
done

mkdir -p mRNA/remove/${PREFIX}
for J in R1 R2;do
  echo "====> $J"
  seqkit grep -j 4 -f <(sed '1d' mRNA/output/${PREFIX}/total_remove_reads_info.tsv | cut -f 1) rRNA_tRNA/filter/${PREFIX}/${J}_conserve_rRNA.fq.gz > mRNA/remove/${PREFIX}/${J}_conserve_rRNA_mRNA.re.fq

  pigz -p 4 mRNA/remove/${PREFIX}/${J}_conserve_rRNA_mRNA.re.fq
done

# ailgn
bash remove_bacteria/script/align.sh remove rRNA_conserve ${PREFIX} rrna
bash remove_bacteria/script/align.sh remove rRNA_tRNA ${PREFIX}  rrna
bash remove_bacteria/script/align.sh remove mRNA ${PREFIX} rrna

bash remove_bacteria/script/align.sh remove rRNA_conserve ${PREFIX} protein_coding
bash remove_bacteria/script/align.sh remove rRNA_tRNA ${PREFIX}  protein_coding
bash remove_bacteria/script/align.sh remove mRNA ${PREFIX} protein_coding

# common
bash remove_bacteria/script/align_rRNA_mRNA.sh remove_align ${PREFIX} rRNA_conserve
bash remove_bacteria/script/align_rRNA_mRNA.sh remove_align ${PREFIX} rRNA_tRNA
bash remove_bacteria/script/align_rRNA_mRNA.sh remove_align ${PREFIX} mRNA
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```

## 6 Visualization
+ plot
```bash
cd remove_bacteria
mkdir -p ../visuliaztion/plot/${PREFIX}

# length
bash script/pre_plot.sh ${PREFIX}
Rscript script/length_hist.r rRNA_conserve ${PREFIX}
Rscript script/length_hist.r rRNA_tRNA ${PREFIX}
Rscript script/length_hist.r mRNA ${PREFIX}

# CG
bash script/pre_plot_CG.sh ${PREFIX}
Rscript script/CG_hist.r rRNA_conserve ${PREFIX}
Rscript script/CG_hist.r rRNA_tRNA ${PREFIX}
Rscript script/CG_hist.r mRNA ${PREFIX}
```

## 7 Make report
+ Raw output
```bash
cd Ath
mkdir -p raw_output/${PREFIX}

# rrna
bowtie2 -p 24 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_rrna \
  -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz \
  -S raw_output/${PREFIX}/rrna.raw.sam \
  2>&1 |
  tee raw_output/${PREFIX}/rrna.bowtie2.log

pigz -p 24 raw_output/${PREFIX}/rrna.raw.sam

# mrna
bowtie2 -p 24 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 -i S,1,0.50 --np 0 \
  --xeq -x index/ath_protein_coding \
  -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz \
  -S raw_output/${PREFIX}/mrna.raw.sam \
  2>&1 |
  tee raw_output/${PREFIX}/mrna.bowtie2.log

pigz -p 24 raw_output/${PREFIX}/mrna.raw.sam

# common
bash remove_bacteria/script/align_rRNA_mRNA.sh raw_output ${PREFIX}
```

```bash
cd Ath

bash remove_bacteria/script/report.sh ${PREFIX}
```

## 8 NJU-seq(rrna)
+ Filter
```bash
PREFIX=Ath_flower_NC

cd Ath
mkdir -p temp/${PREFIX} output/${PREFIX}

pigz -dcf stage_align/${PREFIX}/mRNA/rrna_align.sam.gz |
  parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 4 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
  >temp/${PREFIX}/rrna.out.tmp

parallel -j 4 "
  perl NJU_seq/rrna_analysis/readend_count.pl \
    NJU_seq/data/ath_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 25s 18s 5-8s
# repeat Ath_flower_1 Ath_flower_2 Ath_flower_3
```
+ Score
```bash
TISSUE=flower

parallel -j 4 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/Ath_${TISSUE}_NC/rrna_{}.tsv \\
    output/Ath_${TISSUE}_1/rrna_{}.tsv \\
    output/Ath_${TISSUE}_2/rrna_{}.tsv \\
    output/Ath_${TISSUE}_3/rrna_{}.tsv \\
      >output/Ath_${TISSUE}_rrna_{}_scored.tsv
  " ::: 25s 18s 5-8s
```

## 8 NJU-seq(mrna)
+ Extract reads can't be mapped to rRNA
```bash
PREFIX=Ath_flower_NC

cd Ath
bash NJU_seq/tool/extract_fastq.sh \
  temp/${PREFIX}/rrna.out.tmp \
  data/${PREFIX}/R1.fq.gz data/${PREFIX}/R1.mrna.fq.gz \
  data/${PREFIX}/R2.fq.gz data/${PREFIX}/R2.mrna.fq.gz
```
+ Align
```bash
bowtie2 -p 24 -a -t \
  --end-to-end -D 20 -R 3 \
  -N 0 -L 10 --score-min C,0,0 \
  --xeq -x index/ath_protein_coding \
  -1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz \
  -S stage_align/${PREFIX}/mrna.new.sam \
  2>&1 |
  tee stage_align/${PREFIX}/mrna.new.bowtie2.log

pigz -p 24 stage_align/${PREFIX}/mrna.new.sam
```
+ Filter
```bash
# 双端测序数据匹配了两次，去重，转换为匹配一次
gzip -dcf stage_align/${PREFIX}/mrna.new.sam.gz |
  parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/mrna_analysis/multimatch_judge.pl
  ' | perl NJU_seq/mrna_analysis/multimatch_judge.pl \
  >temp/${PREFIX}/mrna.out.tmp
# bsub -q mpi -n 24 -J "flowerNC" "bash remove_bacteria/script/mrna_filter.sh ${PREFIX}"
head -n 2 temp/${PREFIX}/mrna.out.tmp
# E00517:615:HCJYHCCX2:3:1101:28239:1538  AT5G49540.1     508     12=     ATCTGTTGGGCT
# E00517:615:HCJYHCCX2:3:1101:28239:1538  AT2G17305.1     213     12=     AGCCCAACAGAT

# info
gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
  >data/ath_exon.info
head -n 2 data/ath_exon.info
# 1       3631    3913    +       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT1G01010.1.exon1;rank=1
# 1       3996    4276    +       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=AT1G01010.1.exon2;rank=2

cat temp/${PREFIX}/mrna.out.tmp |
  parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 24 '
    perl NJU_seq/mrna_analysis/dedup.pl \
      --refstr "Parent=transcript:" \
      --transid "AT" \
      --info data/ath_exon.info
  ' |
  perl NJU_seq/mrna_analysis/dedup.pl \
    --refstr "Parent=transcript:" \
    --transid "AT" \
    --info data/ath_exon.info \
    >temp/"${PREFIX}"/mrna.dedup.tmp
# bsub -q mpi -n 24 -J "flowerNC" "bash remove_bacteria/script/dedup.sh ${PREFIX}"

# reads多次匹配,只保留R1匹配情况,reads单次匹配,直接保留
bash NJU_seq/mrna_analysis/almostunique.sh \
  temp/${PREFIX}/mrna.dedup.tmp \
  data/${PREFIX}/R1.mrna.fq.gz \
  temp/${PREFIX} \
  temp/${PREFIX}/mrna.almostunique.tmp
# bsub -q serial -n 2

# 统计基因上不同位点作为起始点和停点出现的次数
perl NJU_seq/mrna_analysis/count.pl \
  temp/${PREFIX}/mrna.almostunique.tmp \
  >temp/${PREFIX}/mrna.count.tmp
# 基因名  停点位置(起始点+长度-2) 匹配reads最后三个碱基  作为起始位点出现次数 作为停点出现次数
# AT5G59400.1     893     G       A       G       0       1
# AT1G21080.3     787     G       A       A       0       1

# 合并位于基因组相同位置的不同转录位点
# 1       3631    3913    +       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT1G01010.1.exon1;rank=1(输入)
gzip -dcf data/ath.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl NJU_seq/mrna_analysis/merge.pl \
    --refstr "Parent=transcript:" \
    --geneid "AT" \
    --transid "AT" \
    -i temp/${PREFIX}/mrna.count.tmp \
    -o output/${PREFIX}/mrna.tsv
# 染色体 位置 正负链 最后三个碱基 基因名称 起始点次数 停点次数
```
+ Calculate valid sequencing depth (average coverage).
```bash
TISSUE=flower
parallel --keep-order -j 4 '
  echo {} >>output/{}/mrna.cov
  bash NJU_seq/presentation/seq_depth.sh \
    temp/{}/mrna.almostunique.tmp \
    output/{}/mrna.tsv \
    >>output/{}/mrna.cov
  ' ::: Ath_${TISSUE}_NC Ath_${TISSUE}_1 Ath_${TISSUE}_2 Ath_${TISSUE}_3
```
+ Score
```bash
perl NJU_seq/mrna_analysis/score.pl \
  output/Ath_${TISSUE}_NC/mrna.tsv \
  output/Ath_${TISSUE}_1/mrna.tsv \
  output/Ath_${TISSUE}_2/mrna.tsv \
  output/Ath_${TISSUE}_3/mrna.tsv \
  >output/Ath_${TISSUE}_mrna_Nm_score.tsv
```

## 8 Statistics and Presentation
+ See the signature of Nm sites
```bash
cd Ath
TISSUE=flower

perl NJU_seq/presentation/signature_count.pl \
  output/Ath_${TISSUE}_mrna_Nm_score.tsv \
  output/Ath_${TISSUE}_mrna_signature.pdf
```

+ End count hist(rrna)
```bash
mkdir visuliaztion/End_count_hist

bash remove_bacteria/script/end_count_hist.sh ${TISSUE} rrna_5-8s
bash remove_bacteria/script/end_count_hist.sh ${TISSUE} rrna_18s
bash remove_bacteria/script/end_count_hist.sh ${TISSUE} rrna_25s
```

+ Ath info
```
1       TAIR10  chromosome      1       30427671        .       .       .       ID=chromosome:1;Alias=Chr1,CP002684.1,NC_003070.9
1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;Name=NAC001;biotype=protein_coding;description=NAC domain containing protein 1 [Source:NCBI gene (formerly Entrezgene)%3BAcc:839580];gene_id=AT1G01010;logic_name=araport11
1       araport11       mRNA    3631    5899    .       +       .       ID=transcript:AT1G01010.1;Parent=gene:AT1G01010;Name=NAC001-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=AT1G01010.1
1       araport11       five_prime_UTR  3631    3759    .       +       .       Parent=transcript:AT1G01010.1
1       araport11       exon    3631    3913    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT1G01010.1.exon1;rank=1
1       araport11       CDS     3760    3913    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    3996    4276    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=AT1G01010.1.exon2;rank=2
1       araport11       CDS     3996    4276    .       +       2       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    4486    4605    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon3;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon3;rank=3
1       araport11       CDS     4486    4605    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    4706    5095    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon4;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon4;rank=4
1       araport11       CDS     4706    5095    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    5174    5326    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon5;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=AT1G01010.1.exon5;rank=5
1       araport11       CDS     5174    5326    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       CDS     5439    5630    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    5439    5899    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon6;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;exon_id=AT1G01010.1.exon6;rank=6
1       araport11       three_prime_UTR 5631    5899    .       +       .       Parent=transcript:AT1G01010.1
```
+ Nm location
```bash
cd Ath
TISSUE=flower
mkdir -p visuliaztion/Nm_location

pigz -dcf data/ath.gff3.gz | grep -v "#" | grep -v "chromosome" > data/ath.gff3
perl remove_bacteria/script/where_is_my_Nm_site.pl data/ath.gff3 output/Ath_${TISSUE}_mrna_Nm_score.tsv | 
  tsv-filter --str-ne 4:mRNA > visuliaztion/Nm_location/${TISSUE}_location.tsv
```
