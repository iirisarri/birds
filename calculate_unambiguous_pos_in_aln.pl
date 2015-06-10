#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;


# Iker Irisarri. University of Konstanz, June 2015
# calculate_unambiguous_pos_in_aln.pl
# modified from random_sample_aln_columns.pl
# Estimates number of unambiguous (a,c,g,t) cols in alignment (fasta format)


my $usage = "calculate_unambiguous_pos_in_aln.pl in.fasta > STDOUT\n";
my $in_fasta = $ARGV[0] or die $usage;


## Store alignment in hash of arrays

my %fasta_hash;
my $header;
my @seq_line;
my @seq;

# using SeqIO for reading in fasta is easier and allows also easily parsing interleaved fasta

my $seqio_obj = Bio::SeqIO->new('-file' => "<$in_fasta", 
								'-format' => "fasta"
							);

while (my $seqio_obj = $seqio_obj->next_seq) {

	my $header = $seqio_obj->primary_id;
	my $seq = $seqio_obj->seq;
	my @seq = split (//, $seq);
	$fasta_hash{$header} = [@seq];
}

#print Dumper \%fasta_hash;


## Get total number of positions in aln & check if seqs are aligned
my $in_aln_total = 0;

foreach my $sp ( keys %fasta_hash ) {
	#print "$sp: ", scalar @{ $fasta_hash{$sp} }, "\n";
	if ( $in_aln_total != 0 && $in_aln_total != scalar @{ $fasta_hash{$sp} } ) {
		print STDERR "WARN: sequence lengths are unequal.
		 They might not be aligned: check your input!\n";
	}
	else {
		$in_aln_total = scalar @{ $fasta_hash{$sp} };
	}
}


## Calculate number of cols consisting of only unambiguous DNA chars

my $defined=0;
my $num_of_seq = scalar keys %fasta_hash;
my $unambiguous_cols=0;

for ( my $i=0; $i<$in_aln_total; $i++ ) {

	# initialize $defined for each aln pos
	$defined = 0;
	
	foreach my $key ( keys %fasta_hash ) {
		#print "$i\t$key\t@{ $fasta_hash{$key} }[$i]\n";
		
		# exit loop if non-DNA char is found
		last if ( @{ $fasta_hash{$key} }[$i] !~ /^(A|C|G|T|a|c|t|g)$/ );
		
		# count if DNA character is found
		if ( @{ $fasta_hash{$key} }[$i] =~ /^(A|C|G|T|a|c|t|g)$/ ) {
			#print "@{ $fasta_hash{$key} }[$i]\n";
			
			$defined++;
		}
		else {
			print STDERR "unexpected error in sequence character!\n";
		}
	}
	if ( $defined == $num_of_seq ) {
	
		$unambiguous_cols++;
	}
}

## Print out info

print "$unambiguous_cols\n";

