#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;


# Iker Irisarri, Feb 2016. University of Konstanz.
# trim_MSA_to_ref_indels.pl
# Trims MSA to positions of reference CDS (zebrafinch)
# This removes 5' and 3' UTR, as well as reading frame shifts if present in sequences (blakcbirds)
# Additionally, checks that codons are complete (divisible by 3) and prints out a warning for files with < 30 codons

my $usage = "trim_MSA_to_ref.pl in.fasta ref_seq > output.fa\n";
my $in_fasta = $ARGV[0] or die $usage;
my $ref_seq = $ARGV[1] or die $usage;

## Store alignment in hash of arrays
open(IN, "<", $in_fasta);

my %fasta_hash;
my $header;
my @seq;
my @pos_to_keep;

=pod
# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$in_fasta",
								'-format' => "fasta");

while ( my $seq_obj = $seqio_obj->next_seq ) {

    $header = $seq_obj->primary_id;
    my $sequence = $seq_obj->seq;
    
    @seq = split (//, $sequence);
    
	$fasta_hash{$header} = [@seq];

}
=cut

# fasta parser without bioperl

while ( my $line =<IN> ) {
	chomp $line;
	if ($line =~ />(.+)/ ) {
		$header = $1;
	}
	else {
		@seq = split (//, $line);
	} 
	$fasta_hash{$header} = [@seq];

}


# get positions without a gap (-) in reference
# These positions will be kept, the rest removed

my @ref_seq = @{ $fasta_hash{$ref_seq} };

my $position = 0;

foreach my $elem ( @ref_seq ) {
	
	$position++;
	
	if ( $elem ne "-" ) {
	
		my $perl_position = $position - 1;
		
		push ( @pos_to_keep, $perl_position );
	}
}

# check that there is at least one position to be kept!
if ( scalar @pos_to_keep == 0 ) {
	print STDERR "WARNING: 0 positions found in reference for $in_fasta\n";
}

my $length = scalar @pos_to_keep;

# print out to stderr if incomplete codons (not divisible by 3)
if ( $length % 3 != 0 ) {
	print STDERR "WARNING: nucleotides not divisible by 3 in $in_fasta\n";
}
# check that there are at least 30 codons (90 nucleotides)
elsif ( $length < 90) {
	print STDERR "WARNING: total number of positions < 90 (< 30 codons)\n";
}
else {
	foreach my $key ( sort keys %fasta_hash ) {

		# print fasta header
		print ">$key\n";

		# print positions to be kept	
		foreach my $p ( @pos_to_keep ) {
	
			print $fasta_hash{$key}[$p];
		}
		print "\n";
	}
}

my $before = scalar @ref_seq;
my $after = scalar @pos_to_keep;

print STDERR "original positions: $before\n";
print STDERR "positions after trimming: $after\n";
