#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

my $usage = "add_id_to_fasta_header.pl fasta id_to_add > new.fasta";
my $fasta = $ARGV[0];
my $info = $ARGV[1];


open(IN, "<", $fasta) or die "Can't open file $fasta!\n";

while ( my $line =<IN> ) {

    chomp $line;
    if ( $line =~ />(.*)/ ) {
	print ">", $info, "_", $1, "\n";
    }
    else {
	print "$line\n";
    }

}

