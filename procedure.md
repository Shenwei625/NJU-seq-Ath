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
    echo "===> $J conserve region align"
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

  perl /script/pre_remove_statistics.pl ../data/${PREFIX}/R1.fq.gz rRNA_conserve/output/${PREFIX}/${J}/${J}_match.tsv \
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
mkdir -p total_rRNA
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
mkdir -p rRNA_tRNA/output/${PREFIX}

find index/rRNA_tRNA -maxdepth 1 -mindepth 1 -type d |
    cut -d "/" -f 3 |
      parallel -j 20 --line-buffer "
        echo '===>{1}'
        mkdir -p rRNA_tRNA/output/${PREFIX}/{1}

        bowtie2 -a -t \
          --end-to-end -D 20 -R 3 \
          -N 0 -L 10 -i S,1,0.50 --np 0 \
          --xeq -x index/rRNA_tRNA/{1}/{1} \
          -1 rRNA_conserve/filter/${PREFIX}/R1_conserve.fq -2 rRNA_conserve/filter/${PREFIX}/R2_conserve.fq \
          -S rRNA_tRNA/output/${PREFIX}/{1}/{1}_align.sam \
          2>&1 |
        tee rRNA_tRNA/output/${PREFIX}/{1}/{1}.bowtie2.log

        pigz rRNA_tRNA/output/${PREFIX}/{1}/{1}_align.sam
      "

# statistics
ls rRNA_tRNA/output/${PREFIX} |
  parallel -j 20 --linebuffer "
  echo '====> {1}'

  pigz -dcf rRNA_tRNA/output/${PREFIX}/{1}/{1}_align.sam.gz |
        grep -v "@" |
        tsv-filter --regex '6:^[0-9]+=$' |
        tsv-filter --str-eq 7:= |
        cut -f 1,10 | sort | uniq > rRNA_tRNA/output/${PREFIX}/{1}/{1}_match.tsv

  perl script/pre_remove_statistics.pl ../data/${PREFIX}/R1.fq.gz rRNA_tRNA/output/${PREFIX}/{1}/{1}_match.tsv \
  rRNA_tRNA/output/${PREFIX}/{1}/reads_info.tsv \
  rRNA_tRNA/output/${PREFIX}/{1}/remove_info.tsv
"

# merge
bash script/merge.sh rRNA_tRNA $PREFIX

tsv-join --filter-file <(cut -d "," -f 1,2 ../../NJU_seq_analysis_ath/bacteria/Bacteria.assembly.collect.csv | tr "," "\t") -H --key-fields name --append-fields Organism_name rRNA_tRNA/output/${PREFIX}/align_statistics.tsv | 
tsv-select -H -f Organism_name --rest last > tem&&
  mv tem rRNA_tRNA/output/${PREFIX}/align_statistics.tsv

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

# align
PREFIX=Ath_flower_NC
mkdir -p mRNA/output/${PREFIX}

find index/mRNA -maxdepth 1 -mindepth 1 -type d |
    cut -d "/" -f 3 |
      parallel -j 20 --line-buffer "
        echo '===>$PREFIX {1}'
        mkdir -p mRNA/output/${PREFIX}/{1}

        bowtie2 -a -t \
          --end-to-end -D 20 -R 3 \
          -N 0 -L 10 -i S,1,0.50 --np 0 \
          --xeq -x index/mRNA/{1}/{1} \
          -1 rRNA_tRNA/filter/${PREFIX}/R1_conserve_rRNA.fq.gz -2 rRNA_tRNA/filter/${PREFIX}/R2_conserve_rRNA.fq.gz \
          -S mRNA/output/${PREFIX}/{1}/{1}_align.sam \
          2>&1 |
        tee mRNA/output/${PREFIX}/{1}/{1}.bowtie2.log

        pigz mRNA/output/${PREFIX}/{1}/{1}_align.sam
      "

# statistics 
bash script/statistics.sh mRNA ${PREFIX}

# merge
bash script/merge.sh mRNA ${PREFIX}

tsv-join --filter-file <(cut -d "," -f 1,2 ../../NJU_seq_analysis_ath/bacteria/Bacteria.assembly.collect.csv | tr "," "\t") -H --key-fields name --append-fields Organism_name rRNA_tRNA/output/${PREFIX}/align_statistics.tsv | 
tsv-select -H -f Organism_name --rest last > tem&&
  mv tem rRNA_tRNA/output/${PREFIX}/align_statistics.tsv

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

