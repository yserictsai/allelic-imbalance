#!/usr/bin/perl

#############################
# 2013.2.3 written by Allen # 
#############################

#die "Usage: perl Parse_mpileup.pl [Combined mpileup results]\n" unless @ARGV; 
#open INFILE, $ARGV[0] or die "Couldn't open infile: $ARGV[0]\n";
#	@combined_mpileup = <INFILE>;
#	chomp (@combined_mpileup);
#close INFILE;


# Tool to report SNP between parents
# usage: 'cat [Combined mpileup results] | ./Parse_mpileup.pl'

print "CDS_name\tSNP_position\tAPO_basecall\tAPO_counts\tIR64_basecall\tIR64_counts\tF1_from_APO\tF1_from_IR64\tF1_total_counts\n";

while ($input = <>) {
        chomp($input);
#for ($i=0; $i<=$#combined_mpileup; $i++) {
	#1: Reference sequence name
	#2: Position
	#3: Reference base
	#4: P1 coverage
	#5: P1 alignment
	#6: P1 mapping quality
	#7: P2 coverage
	#8: P2 alignment
	#9: P2 mapping quality
	#10: F1 coverage
	#11: F1 alignment
	#12: F1 mapping quality

	($ref_name, $ref_pos, $ref_base) = ("", "", "");
	($p1_alignment, $p2_alignment, $f1_alignment) = ("", "", "");
	($p1_A_freq, $p1_C_freq, $p1_G_freq, $p1_T_freq) = (0, 0, 0 ,0);
	($p2_A_freq, $p2_C_freq, $p2_G_freq, $p2_T_freq) = (0, 0 ,0 ,0);
	($f1_A_freq, $f1_C_freq, $f1_G_freq, $f1_T_freq) = (0, 0, 0, 0);
	($p1_basecall, $p2_basecall) = ('N', 'N');

#	($ref_name, $ref_pos, $ref_base, $p1_cov, $p1_alignment, $p2_cov, $p2_alignment, $f1_cov, $f1_alignment) = ($1, $2, $3, $4, $5, $7, $8, $10, $11) if ($combined_mpileup[$i] =~ /(.+)\t(\d+)\t(\w)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)/);	
	($ref_name, $ref_pos, $ref_base, $p1_cov, $p1_alignment, $p2_cov, $p2_alignment, $f1_cov, $f1_alignment) = ($1, $2, $3, $4, $5, $7, $8, $10, $11) if ($input =~ /(.+)\t(\d+)\t(\w)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)/);

	next if($p1_cov == 0 or $p2_cov == 0 or $f1_cov == 0);

	$p1_alignment =~ s/[\.\,]/$ref_base/g;
	@p1_alignment = split('', $p1_alignment);
	for ($j=0; $j<=$#p1_alignment; $j++) {
		if ($p1_alignment[$j] eq '-' or $p1_alignment[$j] eq '+') {
			$j = $j + $p1_alignment[$j+1] + 1;
		} else {
			$p1_A_freq++ if ($p1_alignment[$j] eq 'A');
			$p1_C_freq++ if ($p1_alignment[$j] eq 'C');
			$p1_G_freq++ if ($p1_alignment[$j] eq 'G');
			$p1_T_freq++ if ($p1_alignment[$j] eq 'T');
		}
	}
	$p1_basecall = 'A' if ($p1_A_freq/$p1_cov > 0.8);
	$p1_basecall = 'C' if ($p1_C_freq/$p1_cov > 0.8);
	$p1_basecall = 'G' if ($p1_G_freq/$p1_cov > 0.8);
	$p1_basecall = 'T' if ($p1_T_freq/$p1_cov > 0.8);

        $p2_alignment =~ s/[\.\,]/$ref_base/g;
        @p2_alignment = split('', $p2_alignment);
        for ($j=0; $j<=$#p2_alignment; $j++) {
                if ($p2_alignment[$j] eq '-' or $p2_alignment[$j] eq '+') {
                        $j = $j + $p2_alignment[$j+1] + 1;
                } else {
                        $p2_A_freq++ if ($p2_alignment[$j] eq 'A');
                        $p2_C_freq++ if ($p2_alignment[$j] eq 'C');
                        $p2_G_freq++ if ($p2_alignment[$j] eq 'G');
                        $p2_T_freq++ if ($p2_alignment[$j] eq 'T');
                }
	}
	$p2_basecall = 'A' if ($p2_A_freq/$p2_cov > 0.8);
	$p2_basecall = 'C' if ($p2_C_freq/$p2_cov > 0.8);
	$p2_basecall = 'G' if ($p2_G_freq/$p2_cov > 0.8);
	$p2_basecall = 'T' if ($p2_T_freq/$p2_cov > 0.8);

        $f1_alignment =~ s/[\.\,]/$ref_base/g;
        @f1_alignment = split('', $f1_alignment);
        for ($j=0; $j<=$#f1_alignment; $j++) {
                if ($f1_alignment[$j] eq '-' or $f1_alignment[$j] eq '+') {
                        $j = $j + $f1_alignment[$j+1] + 1;
                } else {
                        $f1_A_freq++ if ($f1_alignment[$j] eq 'A');
                        $f1_C_freq++ if ($f1_alignment[$j] eq 'C');
                        $f1_G_freq++ if ($f1_alignment[$j] eq 'G');
                        $f1_T_freq++ if ($f1_alignment[$j] eq 'T');
                }
	}

	if ($f1_cov >= 5 && $p1_basecall ne 'N' && $p2_basecall ne 'N' && $p1_basecall ne $p2_basecall && $p1_cov >= 5 && $p2_cov >= 5) {
		print "$ref_name\t";
		print "$ref_pos\t";
		print "$p1_basecall\t";
		print "$p1_A_freq\t" if ($p1_basecall eq 'A');
		print "$p1_C_freq\t" if ($p1_basecall eq 'C');
		print "$p1_G_freq\t" if ($p1_basecall eq 'G');
		print "$p1_T_freq\t" if ($p1_basecall eq 'T');
		print "$p2_basecall\t";
		print "$p2_A_freq\t" if ($p2_basecall eq 'A');
		print "$p2_C_freq\t" if ($p2_basecall eq 'C');
		print "$p2_G_freq\t" if ($p2_basecall eq 'G');
		print "$p2_T_freq\t" if ($p2_basecall eq 'T');
		print "$f1_A_freq\t" if ($p1_basecall eq 'A');
		print "$f1_C_freq\t" if ($p1_basecall eq 'C');
		print "$f1_G_freq\t" if ($p1_basecall eq 'G');
		print "$f1_T_freq\t" if ($p1_basecall eq 'T');
		print "$f1_A_freq\t" if ($p2_basecall eq 'A');
		print "$f1_C_freq\t" if ($p2_basecall eq 'C');
		print "$f1_G_freq\t" if ($p2_basecall eq 'G');
		print "$f1_T_freq\t" if ($p2_basecall eq 'T');
		print "$f1_cov\n";
	}
}
