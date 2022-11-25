#! usr/bin/perl
use strict;
use warnings;
use autodie;
use PerlIO::gzip;

=head1 NAME
    pre_remove_statistics.pl  
    -- Statistics before removing reads, including remove Amount(proportion)、length、CG content 
=head1 SYNOPSIS
    perl pre_remove_statistics.pl raw_file.fq.gz remove_reads_name.tsv reads_info.tsv remove_info.tsv
    raw_file.fq.gz          --Original sequencing file
    remove_reads_name.tsv   --File containing the name of reads to be removed, one name per line
    reads_info.tsv          --The name of the output file where the read information is write
    remove_info.tsv         
=cut

open( my $raw_fh1, "<:gzip", $ARGV[0] );
open( my $remove_fh1, "<", $ARGV[1] );
open( my $out_fh, ">", $ARGV[2] );
open( my $out_fh2, ">", $ARGV[3] );

my $dirname;
if ( $ARGV[2] =~ /(.+)\/reads_info.tsv/) {
    $dirname = $1;
}

my @total_reads;
while (<$raw_fh1>) {
    chomp;
    push @total_reads, $_;
    }
close $raw_fh1;

my @remove_name;
while (<$remove_fh1>) {
    chomp;
    push @remove_name, $_;
    }
close $remove_fh1;
my $remove_reads_number = ( $#remove_name + 1 );

my $total_reads_number = 0;
foreach ( @total_reads ) {
    if ( /^@/ ) {
        $total_reads_number++;
    }
}

select $out_fh2;
my $remove_proportion =( ( $remove_reads_number / $total_reads_number ) * 100 );
print "total_reads_number\tmatch_reads_number\tmatch_proportion\n";
print "$total_reads_number\t$remove_reads_number\t";
printf "%.2f%%\n", $remove_proportion;

select $out_fh;
# reads information
my $i = 0;
my @discard_reads;
print "reads_name\treads_sequence\treads_length\tCG_content\n";
foreach ( @remove_name ) {
    my $reads_name;
    my $reads_seq;
    if ( /^(\S+)\t(\S+)/ ) {
        $reads_name = $1;
        $reads_seq = $2;
    }
    my $reads_length = length($reads_seq);

    my $C_count = $reads_seq =~ tr/C/C/;
    my $G_count = $reads_seq =~ tr/G/G/;
    my $CG_count = ( ( $C_count + $G_count ) / $reads_length );

    print "$reads_name\t$reads_seq\t$reads_length\t$CG_count\n";
}


__END__
