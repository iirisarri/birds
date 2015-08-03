#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

# Iker Irisarri, Aug 2015. University of Konstanz
# Modif from extract_from_fasta_by_name.pl
# script reads tab-delimited file of RBBH pairs (generated with parse_blast_RBBH.pl),
# finds file named with name of trinity component and appends there the sequence of the RBBH

my $usage = "rbbh_2_fasta_file.pl source_fasta query_file > STDOUT\n";
my $fasta = $ARGV[0] or die $usage; #"Taeniopygia_guttata.taeGut3.2.4.cds.all.tid.fa";
my $query = $ARGV[1] or die $usage;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
								'-format' => "fasta");
        
# declare hash
my %queries;
my %in_seqs;

print STDERR "reading in RBBH from file...\n";

#open query file, chomp line and save it into the array with push
#the hash will contain gene names (variable $line) as keys and 1 as a value (that's random)
open (IN , $query) or die "Can't open $query, $!\n";

while (my $line = <IN>){
    chomp $line;

    my ($comp, $enst) = split ( "\t", $line );

    $queries{$comp} = $enst;

}

print STDERR "reading fasta...\n";

while (my $seq_obj = $seqio_obj->next_seq) {

    my $seqname = $seq_obj->primary_id;
    my $sequence = $seq_obj->seq;
    
    $in_seqs{$seqname} = $sequence;

}

print STDERR "finding correspondence and printing to files...\n";

foreach my $key ( keys %queries ) {

    my $zebrafinch = $queries{$key};

    if ( exists $in_seqs{$zebrafinch} ) {

	my $file = $key . "_consensus.fa";

	open (OUT, ">>", $file);

	print "printing $zebrafinch to $file!\n";

	print OUT  ">zebrafinch\n";
	print OUT "$in_seqs{$zebrafinch}\n";
    }

}

