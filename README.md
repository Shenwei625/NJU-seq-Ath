# NJU-seq in *Arabidopsis thaliana*

## Reference and index
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

## Data Selection and quality overview
```bash
# Select data for analysing
ID='NJU6220'
PREFIX='Ath_root_RF_NC'

mkdir -p "data/${PREFIX}" "temp/${PREFIX}" "output/${PREFIX}"

# cut adapter


# Quality control for clean data
# Turn pair-end sequence data to single-end file
perl NJU_seq/quality_control/pe_consistency.pl \
  data/"${PREFIX}"/NJU6220_S1_L002_R1_001.fastq.gz data/"${PREFIX}"/NJU6220_S1_L002_R2_001.fastq.gz \
  temp/"${PREFIX}".fastq.gz
```



