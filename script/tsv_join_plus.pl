#! usr/bin/perl
use strict;
use warnings;
use autodie;

=head1 NAME
    tsv_join_plus.pl
    -- According to the first columns of two files, add another column of file2 to file1, Lines not saved in file2 will be replaced with 0
=head1 SYNOPSIS
    tsv_join_plus.pl raw_file append_file
=cut

open( my $in_fh1, "<", $ARGV[0] );
open( my $in_fh2, "<", $ARGV[1] );

my %info1;
while (<$in_fh1>) {
    if ( ! (/^reads_length|^CG_content/)) {
        chomp;
        my ( $dkey, $dvalue ) = split(/\t/, $_);
        $info1{$dkey} = $dvalue;
    }
}
close $in_fh1;

my %info2;
while (<$in_fh2>) {
    if ( ! (/^reads_length|^CG_content/)) {
        chomp;
        my ( $dkey, $dvalue ) = split(/\t/, $_);
        $info2{$dkey} = $dvalue;
    }
}
close $in_fh2;

my @filekey;
foreach ( keys %info2 ) {
    push @filekey, $_;
}

print "reads_length\tTotal_number\tKeep_number\tKeep_proportion\n";
foreach my $key ( sort { $a <=> $b } keys %info1 ) {
    if ( grep { $key eq $_ } @filekey ) {
        my $remove_proportion = ( ( $info2{$key} / $info1{$key} ) * 100 );
        my $rm_print =`printf "%.2f%%" $remove_proportion`;
        print "$key\t$info1{$key}\t$info2{$key}\t$rm_print\n";
    }else {
        print "$key\t$info1{$key}\t0\t0\n";
    }
}
__END__
