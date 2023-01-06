#!/usr/bin/perl

use strict;
use warnings;

die "Usage: perl Parse_parents_BAM.pl [BAM] \n" unless @ARGV; 

my %p1_pos2base = ();
my %p2_pos2base = ();
my $bam_file = "$ARGV[0]";
my @sam = ();
my %loc_2_length = ();

open (P_SNP, 'all_rep_snpcall.txt');
while(<P_SNP>) {
	my ($loc, $ref_pos, $p1_base, $p2_base) = ($1, $2, $3, $4) if ($_ =~ /(.+?)\t(\d+)\t\w\t(\w)\t(\w)\W?/);  #why .+? \W?
	$p1_pos2base{$loc}[$ref_pos] = "$p1_base";
	$p2_pos2base{$loc}[$ref_pos] = "$p2_base";
}
close P_SNP;


open ( P1, ">Rep1_from_apo_c.sam" );
open ( P2, ">Rep1_from_ir64_c.sam" );
open ( CN, ">Rep1_canot_distin_c.sam" );
#open ( NOSNP, ">nosnp.sam" );





#open BAM, "samtools view -h $bam_file | sort --temporary-directory=`pwd` -k 3,3 |";
open BAM, "samtools view -h $bam_file |";
my ($tot_match, $tot_insertion, $tot_deletion, $tot_soft_clipping, $tot_pairs, $tot_tlen) = (0, 0, 0, 0, 0, 0);
while(<BAM>) {
#	next if(/^(\@)/);	## Skipping the header lines (if you used -h in the samools command)
	if (/^\@/) {	## If headers
		$loc_2_length{$1} = "$2" if ($_ =~ /\@SQ\tSN:(.+?)\tLN:(\d+)/);      ###why ?
	} else {
		s/\n//;  s/\r//;	## Removing new line
		@sam = split(/\t+/);	## Splitting SAM line into array
		my ($qname, $flag, $rname, $pos, $cigar, $rnext, $pnext, $tlen, $seq) = ($sam[0], $sam[1], $sam[2], $sam[3], $sam[5], $sam[6], $sam[7], $sam[8], $sam[9]);
		my ($match, $insertion, $deletion, $soft_clipping, $tot_read_length) = (0, 0, 0, 0, 0);

		next if ($rname =~ /\*/);	## Ignore unmapped reads

		## From the CIGAR string, fill in the insertions and deletions
		my $cigar_ori = $cigar;
		while ($cigar !~ /^$/){
			if ($cigar =~ /^([0-9]+[MIDS])/){
				my $cigar_part = $1;
				if ($cigar_part =~ /(\d+)M/){
					$match += $1;
				} elsif ($cigar_part =~ /(\d+)I/){
					$insertion += $1;
				} elsif ($cigar_part =~ /(\d+)D/){
					$deletion += $1;
				} elsif ($cigar_part =~ /(\d+)S/){
					$soft_clipping += $1;
				}
				$cigar =~ s/$cigar_part//;     # remove cigar_part?
			} elsif ($cigar =~ /\*/){
				last;
			} else {
				die "Unexpected cigar: $qname $cigar\n";
			}
		}
#		$tot_read_length = $match + $insertion + $deletion + $soft_clipping;
		$tot_read_length = $match;

		my ($f1_same_with_p1, $f1_same_with_p2) = (0, 0);
		my (@f1_haplotype_ary, @p1_haplotype_ary, @p2_haplotype_ary) = ();
		for (my $read_pos=0; $read_pos<$tot_read_length; $read_pos++) {
			my $read_on_ref_pos = $read_pos + $pos;    ##first bp to ref
			if (defined(${$p1_pos2base{$rname}}[$read_on_ref_pos]) && defined(${$p2_pos2base{$rname}}[$read_on_ref_pos])) {       # in snpcall's base
				my $p1_base = ${$p1_pos2base{$rname}}[$read_on_ref_pos];
				my $p2_base = ${$p2_pos2base{$rname}}[$read_on_ref_pos];
				my $f1_base = substr($seq, $read_pos, 1);
#				print "$qname, $read_on_ref_pos, F1:$f1_base, P1:$p1_base, P2:$p2_base\n";
				$f1_same_with_p1++ if ($f1_base eq $p1_base);
				$f1_same_with_p2++ if ($f1_base eq $p2_base);
				push (@f1_haplotype_ary, $f1_base);
				push (@p1_haplotype_ary, $p1_base);
				push (@p2_haplotype_ary, $p2_base);
			} #else {
			#	print NOSNP "$_\n";
			#}
		}
		my $f1_haplotype = join('', @f1_haplotype_ary);
		my $p1_haplotype = join('', @p1_haplotype_ary);
		my $p2_haplotype = join('', @p2_haplotype_ary);
		if ($f1_same_with_p1 > 0 && $f1_haplotype eq $p1_haplotype) {
			#print  P1 "$qname\tRead from P2\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
			print  P1 "$_\n";
		} elsif ($f1_same_with_p2 > 0 && $f1_haplotype eq $p2_haplotype) {	
			#print  P2 "$qname\tRead from P2\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
			print  P2 "$_\n";
		} elsif ($f1_haplotype ne $p1_haplotype && $f1_haplotype ne $p2_haplotype) {
			#print  CN "$qname\tRead cannot distinct\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
 			print  CN "$_\n";
		}
		($match, $insertion, $deletion, $soft_clipping) = (0, 0, 0, 0);
	}
}
close BAM;
close P1;
close P2;
close CN;
#close NOSNP;


sub round {
	my $val = shift;
	my $col = shift;
	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;
	return int($val * $r + $a) / $r;
}
