#!/usr/bin/perl

use strict;
use warnings;

die "Usage: perl Parse_parents_BAM.pl [BAM] \n" unless @ARGV; 

my ($bam_file, $snpcall_file) = @ARGV;
my %p1_pos2base = ();
my %p2_pos2base = ();

open P_SNP, "<$snpcall_file";
while(<P_SNP>) {
	my ($loc, $ref_pos, $p1_base, $p2_base) = ($1, $2, $3, $4) if ($_ =~ /(.+)\t(\d+)\t\w\t(\w)\t(\w)/);  
	$p1_pos2base{$loc}[$ref_pos] = "$p1_base";
	$p2_pos2base{$loc}[$ref_pos] = "$p2_base";
}
close P_SNP;


#print $p1_pos2base{"LOC_Os01g01920.2"}[1] # it's a hash array which index = refpose has a base other index is empty 
open ( P1, ">Rep1_from_apo_s.sam" );
open ( P2, ">Rep1_from_ir64_s.sam" );
open ( CN, ">Rep1_canot_distin_s.sam" );
#open ( NOSNP, ">nosnp.sam" );
#open BAM, "samtools view -h $bam_file | sort --temporary-directory=`pwd` -k 3,3 |";
open BAM, "samtools view -h $bam_file |";

my @cigar_array = ();
my @sam = ();

while(<BAM>) {
	chomp;
        my $line = $_;
	if( $line =~ /^(\@)/){
		next;	## Skipping the header lines (if you used -h in the samools command)
	}else {
		s/\n//;  s/\r//;	## Removing new line
		@sam = split(/\t+/);	## Splitting SAM line into array
		my ($qname, $flag, $rname, $pos, $cigar, $rnext, $pnext, $tlen, $seq) = ($sam[0], $sam[1], $sam[2], $sam[3], $sam[5], $sam[6], $sam[7], $sam[8], $sam[9]);
		my ($match, $insertion, $deletion, $tot_length) = (0, 0, 0, 80);
		next if ($rname =~ /\*/);	## Ignore unmapped reads
		my ($f1_same_with_p1, $f1_same_with_p2) = (0, 0);
                my (@f1_haplotype_ary, @p1_haplotype_ary, @p2_haplotype_ary) = ();		
		my $index = 0;
		my $fl = 0;
		my $pos_index = 0;
		@cigar_array = ( $cigar =~ /(\d+[MID])/g );

		for ( my $i= 0; $i<=$#cigar_array; $i++ ){
			my $insert_f1 = 0;
			my $deletion_fl = 0;
			$index = $fl;
			if ($cigar_array[$i] =~ /(\d+)M/){
				$match += $1;
			} elsif ($cigar_array[$i] =~ /(\d+)I/){
				$insertion += $1;
				$insert_f1 = 1;
			} elsif ($cigar_array[$i] =~ /(\d+)D/){
				$deletion += $1;
				$deletion_fl = 1;
			}
			$fl = $match + $insertion;
			$pos_index = $deletion - $insertion;
	
			for ( my $read_pos=0+$index; $read_pos<$tot_length; $read_pos++){
				last if ( $insert_f1 == 1 or $deletion_fl == 1 );
				my $read_on_ref_pos = $read_pos + $pos;
				
				if (defined(${$p1_pos2base{$rname}}[$read_on_ref_pos+$pos_index]) && defined(${$p2_pos2base{$rname}}[$read_on_ref_pos+$pos_index])) {            
					my $p1_base = ${$p1_pos2base{$rname}}[$read_on_ref_pos+$pos_index];
					my $p2_base = ${$p2_pos2base{$rname}}[$read_on_ref_pos+$pos_index];
					my $f1_base = substr($seq, $read_pos, 1);
#                               print "$qname, $read_on_ref_pos, F1:$f1_base, P1:$p1_base, P2:$p2_base\n";
					$f1_same_with_p1++ if ($f1_base eq $p1_base);
					$f1_same_with_p2++ if ($f1_base eq $p2_base);
					push (@f1_haplotype_ary, $f1_base);
					push (@p1_haplotype_ary, $p1_base);
					push (@p2_haplotype_ary, $p2_base);
				}			
				last if ( $read_pos+1 == $fl );
			}
		}
		my $f1_haplotype = join('', @f1_haplotype_ary);
		my $p1_haplotype = join('', @p1_haplotype_ary);
		my $p2_haplotype = join('', @p2_haplotype_ary);
		if ($f1_same_with_p1 > 0 && $f1_haplotype eq $p1_haplotype) {
			#print "$qname\t$rname\t$pos\t$cigar\tRead from P1\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
			print  P1 "$_\n";
		} elsif ($f1_same_with_p2 > 0 && $f1_haplotype eq $p2_haplotype) {	
			#print "$qname\t$rname\t$pos\t$cigar\tRead from P2\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
			print  P2 "$_\n";
		} elsif ($f1_haplotype ne $p1_haplotype && $f1_haplotype ne $p2_haplotype) {
			#print  CN "$qname\tRead cannot distinct\tF1 Haplotype: $f1_haplotype, P1 Haplotype: $p1_haplotype, P2 Haplotype: $p2_haplotype\tP1: $f1_same_with_p1, P2: $f1_same_with_p2\n";
 			print  CN "$_\n";
		}
		#($match, $insertion, $deletion, $soft_clipping) = (0, 0, 0, 0);
	}
}



#close BAM;
#close P1;
#close P2;
#close CN;
#close NOSNP;




