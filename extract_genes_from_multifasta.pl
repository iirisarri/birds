#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;
use Data::Dumper;

my $usage = "extractSeqFromFasta fasta_file > stdout\n";
my $fasta = $ARGV[0] or die $usage;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
				'-format' => "fasta");

my %fasta_hash;                   

while (my $seq_obj = $seqio_obj->next_seq){

    $seq_obj->primary_id =~ /(B\d\d)_(.*)/;
    my $individual = $1;
    my $transcript = $2;
    my $seq = $seq_obj->seq;

    $fasta_hash{$transcript}{$individual} = $seq;

}

#print STDERR Dumper \%fasta_hash;

foreach my $key ( keys %fasta_hash ) {

    my $outfile = $key . ".fa";
    open (OUT, ">", $outfile);

    my %inner_hash = %{ $fasta_hash{$key} };

    foreach my $indiv ( sort keys %inner_hash ) {
	print OUT ">", $indiv, "_", $key, "\n";
	print OUT "$inner_hash{$indiv}\n";
    }
    close(OUT);
}

