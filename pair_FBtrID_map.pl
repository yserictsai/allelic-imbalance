#!/usr/bin/perl -w

use strict;

die "

	Example: $0 correlation.without_less_than_1.txt_cut_0.95_pairs.txt /CSBBNAS/Project/NGS_special/Manuscript2/mappingtable

" if !@ARGV;

my ($file, $map) = @ARGV;
my %mapping = ();
open FI, "$map";
my $head = <FI>;
while (<FI>){
	chomp $_;
	my @t = split ("\t", $_);
	$mapping{uc($t[$#t])} = $t[0];
	$mapping{uc($t[($#t-1)])} = $t[0];
}
close FI;

my %pair = ();
open FI, "$file";
while (<FI>){
	my @t = split(/\s+/, $_);
	if (defined $mapping{uc($t[0])} && defined $mapping{uc($t[1])}){
		my @tmp = ($mapping{uc($t[0])}, $mapping{uc($t[1])});
		@tmp = sort @tmp;
		my $line = join(" ", @tmp);
		$pair{$line} = $t[2];
	}
}
close FI;

foreach my $e (sort keys %pair){
	print "$e\t$pair{$e}\n"
}
