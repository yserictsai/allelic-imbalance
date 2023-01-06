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
open ( P1, ">mod_Rep2_extract_APO_s.sam" );
#open ( P2, ">tt_Rep1_from_ir64_c.sam" );
#open ( CN, ">tt_Rep1_canot_distin_c.sam" );
#open ( NOSNP, ">nosnp.sam" );
#open BAM, "samtools view -h $bam_file | sort --temporary-directory=`pwd` -k 3,3 |";
open BAM, "samtools view -h $bam_file |";

my @cigar_array = ();
my @sam = ();
my $k = 0;
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
		my $break_flag = 0;
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
				$break_flag = 1 if (defined(${$p1_pos2base{$rname}}[$read_on_ref_pos+$pos_index]));
				last if ( $read_pos+1 == $fl or $break_flag == 1 );
			}
			last if ( $break_flag == 1 );
		}
		print P1 "$_\n" if ( $break_flag == 1 );
#		print P2 "$_\n" if ( $break_flag == 0 );

	}
	#$k = $k + 1 ;


}

#print $k;

close BAM;
close P1;
#close P2;
#close CN;
#close NOSNP;




