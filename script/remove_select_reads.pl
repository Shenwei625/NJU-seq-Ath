#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;

=head1 NAME
remove_select_reads.pl  -- Remove selected reads 
=head1 SYNOPSIS
    perl remove_select_reads.pl -r file1.fastq -f file2.tsv
        options:
            --help\-h   Help message
            --raw\-r    The raw fastq file(compressed) that you want to process
            --file\-f   File containing reads name that you want to remove(One name per line, No prefix @ before name) 
            --stdin     Get FastA from STDIN 
=cut

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'raw|r=s'   => \my $raw_file,
    'file|f=s'  => \my $reads_names,
    'stdin'     => \my $stdin,
) or Getopt::Long::HelpMessage(1);

my @ALL;
if ( defined($raw_file) ) {
    open( my $in_fh1, "<", $raw_file );
    while (<$in_fh1>) {
        chomp;
        push @ALL, $_;
    }
    close $in_fh1;
}elsif ( defined($stdin) ) {
    @ALL = <>;
}else {
    die("You should provide raw data!")
}
chomp(@ALL);

my @total_names;
if ( defined($reads_names) ) {
    open( my $FA, "<", $reads_names );
    while (<$FA>) {
        chomp;
        push @total_names, $_;
    }
    close $FA;
}
else {
    die("You should provide reads names file!")
}

my $i = 0;
while ( $i < $#ALL ) {
    if ( $ALL[$i] =~ /^@(\S+)/ ) {
        my $name_index = $i;
        my $f_name = $1;
        $i++;
        my $seq1 = $ALL[$i];
        $i++;
        my $info = $ALL[$i];
        $i++;
        my $quality = $ALL[$i];
        $i++;
        print "$ALL[$name_index]\n$seq1\n$info\n$quality\n" if( ! ( "@total_names" =~ m/$f_name/ ));
    }else {
        exit;
    }
}

__END__