#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub COUNT_CIGAR {
    my $CIGAR   = shift;
    my $LENGTH  = 0;
    my $MATCHED = 0;
    my $INSERT  = 0;
    my $DELETE  = 0;
    my $XMATCH  = 0;

    # Alignment column containing two identical letters
    while ( $CIGAR =~ /([0-9]+)=/g ) {
        $LENGTH  += $1;
        $MATCHED += $1;
    }

    # Alignment column containing a mismatch, i.e. two different letters
    while ( $CIGAR =~ /([0-9]+)X/g ) {
        $LENGTH += $1;
        $XMATCH += $1;
    }

    # Deletion (gap in the target sequence)
    while ( $CIGAR =~ /([0-9]+)D/g ) {
        $LENGTH += $1;
        $DELETE += $1;
    }

    # Insertion (gap in the query sequence)
    while ( $CIGAR =~ /([0-9]+)I/g ) {
        $INSERT += $1;
    }

    return ( $LENGTH, $MATCHED, $XMATCH, $DELETE, $INSERT );
}

while (<>) {
    my ( $read_name, $cigar ) = split /\t/;
    my ( $l, $m, $x, $d, $i ) = COUNT_CIGAR($cigar);

    if ( $l == $m ) {
        print "$read_name\n";
    }
}

