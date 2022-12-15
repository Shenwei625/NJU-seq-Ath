#! usr/bin/bash
# bash pre_plot.sh PREFIX
$PREFIX=$1

J=../data/${PREFIX}
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

for D in rRNA_conserve rRNA_tRNA mRNA;do
  sed '1d' $dir/output/${PREFIX}/total_remove_reads_info.tsv >> ${PREFIX}_remove.tsv
done

cut -f 3 ${PREFIX}_remove.tsv | 
  sed '1d' |
  sort -n | uniq -c |
  perl -ne'
    s/^\s+//;
    print "$_";
  ' |
  tr " " "\t" | tsv-select --fields 2,1 > ${PREFIX}_hist.tsv

(echo -e "reads_length\tRemove_number" && cat ${PREFIX}_hist.tsv) > tem&&
  mv tem ${PREFIX}_hist.tsv

perl script/tsv_join_plus.pl ../data/${PREFIX}/length_distribution.tsv ${PREFIX}_hist.tsv > tem&&
  mv tem ../data/${PREFIX}/length_distribution.tsv

echo -e "reads_length\tnumber\tgroup" > ../data/${PREFIX}/plot.tsv
sed '1d' ../data/${PREFIX}/length_distribution.tsv | perl -ne'
  chomp;
  if (/^(\S+)\t(\S+)\t(\S+)/) {
    my $reads_length = $1;
    my $total_number = $2;
    my $remove_number = $3;
    my $keep_number = ( $total_number - $remove_number );
    print "$reads_length\t$keep_number\tKeep\n";
    print "$reads_length\t$remove_number\tRemove\n";
  }
' >> ../data/${PREFIX}/plot.tsv

