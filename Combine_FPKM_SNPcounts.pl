#!/usr/bin/perl

# Combine_FPKM_SNPcounts.pl
# AUTHOR: Allen Kao
# LAST REVISED: 2013.2.4
# CONTACT: kaosm16@gmail.com

# Tool to combine FPKM and SNP calls
# usage: 'perl Combine_FPKM_SNPcounts.pl [Combined FPKM table] [Parsed mpileup result]'

die "Usage: perl Combine_FPKM_SNPcounts.pl [Combined FPKM table] [Parsed mpileup result]\n" unless @ARGV; 

# Open FPKM table
open INFILE, $ARGV[0] or die "Couldn't open infile: $ARGV[0]\n";
	my @combined_fpkm = <INFILE>;
	chomp (@combined_fpkm);
close INFILE;
# Open SNP counts table
open INFILE, $ARGV[1] or die "Couldn't open infile: $ARGV[1]\n";
        my @combined_mpileup = <INFILE>;
        chomp (@combined_mpileup);
close INFILE;

# Declare FPKM and counts hash
my (%p1_counts_hash, %p1_fpkm_hash);
my (%p2_counts_hash, %p2_fpkm_hash);
my (%f1_counts_hash, %f1_fpkm_hash);

# Make a FPKM and counts hash (keys=CDS, value=FPKM, counts)
for (my $i=0; $i<=$#combined_fpkm; $i++) {
	next if($combined_fpkm[$i] =~ /^CDS_name/);
	my ($cds, $p1_counts, $p1_fpkm, $p2_counts, $p2_fpkm, $f1_counts, $f1_fpkm) = ($1, $2, $3, $4, $5, $6, $7) if ($combined_fpkm[$i] =~ /(.+)\t(\d+\.\d+)\t(\d+\.\d+)\t(\d+\.\d+)\t(\d+\.\d+)\t(\d+\.\d+)\t(\d+\.\d+)/);
	$p1_counts_hash{$cds} = round($p1_counts, 0);
	$p1_fpkm_hash{$cds} = "$p1_fpkm";
	$p2_counts_hash{$cds} = round($p2_counts, 0);
	$p2_fpkm_hash{$cds} = "$p2_fpkm";
	$f1_counts_hash{$cds} = round($f1_counts, 0);
	$f1_fpkm_hash{$cds} = "$f1_fpkm";
}

# Declare fields in SNP counts table
my ($cds_name, $snp_pos) = ("", 0);
my ($p1_basecall, $p1_snp_counts) = ("", 0);
my ($p2_basecall, $p2_snp_counts) = ("", 0);
my ($f1_from_p1, $f1_from_p2, $f1_from_parents, $f1_tot_counts) = (0, 0, 0, 0);

# Print header
print "\"CDS\",";
print "\"SNP_position\",";
print "\"APO_CDS_FPKM\",";
print "\"APO_CDS_counts\",";
print "\"APO_SNP_basecall\",";
print "\"APO_SNP_counts\",";
print "\"IR64_CDS_FPKM\",";
print "\"IR64_CDS_counts\",";
print "\"IR64_SNP_basecall\",";
print "\"IR64_SNP_counts\",";
print "\"F1_CDS_FPKM\",";
print "\"F1_CDS_counts\",";
print "\"F1_from_APO\",";
print "\"F1_from_IR64\",";
print "\"F1_total_SNP_counts\"\n";

# Add CDS counts and FPKM to each SNP call
for (my $i=0; $i<=$#combined_mpileup; $i++) {
	#1: CDS name
	#2: SNP position
	#3: P1 basecall
	#4: P1 SNP counts
	#5: P2 basecall
	#6: P2 SNP counts
	#7: F1 from P1
	#8: F1 from P2
	#9: F1 total SNP counts

	next if($combined_mpileup[$i] =~ /^CDS_name/);
	($cds_name, $snp_pos, $p1_basecall, $p1_snp_counts, $p2_basecall, $p2_snp_counts, $f1_from_p1, $f1_from_p2, $f1_tot_counts) = ($1, $2, $3, $4, $5, $6, $7, $8, $9) if ($combined_mpileup[$i] =~ /^(.+)\t(\d+)\t([ACGT])\t(\d+)\t([ACGT])\t(\d+)\t(\d+)\t(\d+)\t(\d+)$/);	

	$f1_from_parents = $f1_from_p1 + $f1_from_p2;
	if ($f1_from_parents <= $f1_tot_counts) {
		print "\"$cds_name\",";
		print "\"$snp_pos\",";
		print "\"$p1_fpkm_hash{$cds_name}\",";
		print "\"$p1_counts_hash{$cds_name}\",";
		print "\"$p1_basecall\",";
		print "\"$p1_snp_counts\",";
		print "\"$p2_fpkm_hash{$cds_name}\",";
		print "\"$p2_counts_hash{$cds_name}\",";	
		print "\"$p2_basecall\",";
		print "\"$p2_snp_counts\",";
		print "\"$f1_fpkm_hash{$cds_name}\",";
		print "\"$f1_counts_hash{$cds_name}\",";
		print "\"$f1_from_p1\",";
		print "\"$f1_from_p2\",";
		print "\"$f1_tot_counts\"\n";
	}
}


sub round {
        my $val = shift;
        my $col = shift;
        my $r = 10 ** $col;
        my $a = ($val > 0) ? 0.5 : -0.5;
        return int($val * $r + $a) / $r;
}
