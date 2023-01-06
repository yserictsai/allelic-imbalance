#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

die "Usage: perl Temp.pl [Contigs FASTA]\n" unless @ARGV;

my ($ref_file, $basecall) = @ARGV;
my %ref_seq = ();



open IN2, "<$basecall";


### Read sequence in FASTA of both reference and contigs   to be a hash key = id value = seq
my $seqio_obj = Bio::SeqIO -> new(-file => "$ref_file", -format => "fasta");
while (my $seq_obj = $seqio_obj -> next_seq) {
        my $id = $seq_obj -> display_id;
        my $seq = $seq_obj -> seq;
        die "The contigs contain irregular symbols in $id!\n" if ($seq =~ /\W/);
        $ref_seq{$id} = "$seq";
}

my @ref_name;
my @pos;
my @indel;

while( <IN2> ){
	chomp;
	my $line = $_;
	if ( $line =~ /(.+)\t(\d+)\t(\S+)/){
		push(@ref_name, $1);
		push(@pos, $2);
		push(@indel, $3);
	}
}
close (IN2);


my $i;
my $sign;
my $num;
my $seq;
for ( $i=0; $i <=$#ref_name; $i++ ){
	if ( $indel[$i] =~  /([+-])(\d+)([A-Z]+)/ ){
		$sign = $1;
		$num = $2;
		$seq = $3;
	}
	substr($ref_seq{$ref_name[$i]}, $pos[$i], 0) = $seq if ( $sign eq '+');
	substr($ref_seq{$ref_name[$i]}, $pos[$i], $num) = '' if ( $sign eq '-');
}



for my $index (keys %ref_seq){
	print ">$index\n";
	print_sequence($ref_seq{$index}, 60);
}



sub print_sequence {
    my($sequence, $length) = @_;
    use strict;
    use warnings;
    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
        print substr($sequence, $pos, $length), "\n";
	}
}
=cut
