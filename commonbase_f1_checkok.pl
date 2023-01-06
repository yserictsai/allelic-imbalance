#!/usr/bin/perl


#use strict;
#use warnings;

#print "CDS_name\tSNP_position\tAPO_C_basecall\tAPO_S_basecall\tAPO_C_counts\tAPO_S_counts\tIR64_C_basecall\tIR64_S_basecall\tIR64_C_counts\tIR64_S_counts";


my $input;
my $ref_name;
my $ref_pos;
my $ref_base;
my $p1a_cov;
my $p1a_alignment;
my $p1b_cov;
my $p1b_alignment;
#my $p2a_cov;
#my $p2a_alignment;
#my $p2b_cov;
#my $p2b_alignment;
my $p1a_basecall;
my $p1b_basecall;
#my $p2a_basecall;
#my $p2b_basecall;
my @p1a_alignment;
my @p1b_alignment;
#my @p2a_alignment;
#my @p2b_alignment;
my $p1a_A_freq;
my $p1a_C_freq;
my $p1a_G_freq;
my $p1a_T_freq;
my $p1b_A_freq;
my $p1b_C_freq;
my $p1b_G_freq;
my $p1b_T_freq;
#my $p2a_A_freq;
#my $p2a_C_freq;
#my $p2a_G_freq;
#my $p2a_T_freq;
#my $p2b_A_freq;
#my $p2b_C_freq;
#my $p2b_G_freq;
#my $p2b_T_freq;
my $p1_mixbasecall;
#my $p2_mixbasecall;


while ( $input = <> ) {
        chomp($input);
        $p1a_basecall = 'X';
        $p1b_basecall = 'X';
#        $p2a_basecall = 'X';
 #       $p2b_basecall = 'X';
        my $index1 = 0;
        my $index2 = 0;
        if ($input =~ /(.+)\t(\d+)\t(\w)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)/){
                $ref_name = $1;
                $ref_pos  = $2;
                $ref_base = $3;
                $p1a_cov  = $4;
                $p1a_alignment = $5;
                $p1b_cov = $7;
                $p1b_alignment = $8;
        }
        # print "$ref_name\n$ref_pos\n$ref_base\n$p1a_cov\n$p1a_alignment\n$p1b_cov\n$p1b_alignment\n$p2a_cov\n$p2a_alignment\n$p2b_cov\n$p2b_alignment\n";


### apo_c_basecall ###


        $p1a_alignment =~ s/[\.\,]/$ref_base/g;
        @p1a_alignment = split('', $p1a_alignment);
        #print "@p1a_alignment\n";
        #print "$p1a_cov\n";
        $p1a_A_freq = 0;
        $p1a_C_freq = 0;
        $p1a_G_freq = 0;
        $p1a_T_freq = 0;
        if( $p1a_cov >= 5 ){
                for ( my $j=0; $j<=$#p1a_alignment; $j++  ) {
			if ( $p1a_alignment[$j] eq '-' or $p1a_alignment[$j] eq '+' ){
				$j = $j + $p1a_alignment[$j+1] + 2;
        		} else {        
		        $p1a_A_freq++ if ($p1a_alignment[$j] =~ /[Aa]/ );
                        $p1a_C_freq++ if ($p1a_alignment[$j] =~ /[Cc]/ );
                        $p1a_G_freq++ if ($p1a_alignment[$j] =~ /[Gg]/ );
                        $p1a_T_freq++ if ($p1a_alignment[$j] =~ /[Tt]/ );
	                }
		}
                #print "$p1a_A_freq\t$p1a_C_freq\t$p1a_G_freq\t$p1a_T_freq\n";
                $p1a_basecall = 'A' if ($p1a_A_freq/$p1a_cov > 0.9);
                $p1a_basecall = 'C' if ($p1a_C_freq/$p1a_cov > 0.9);
                $p1a_basecall = 'G' if ($p1a_G_freq/$p1a_cov > 0.9);
                $p1a_basecall = 'T' if ($p1a_T_freq/$p1a_cov > 0.9);

        }
        #print "$ref_name\t$ref_pos\t$p1a_basecall\n";


### apo_s_basecall ###
        $p1b_alignment =~ s/[\.\,]/$ref_base/g;
        @p1b_alignment = split('', $p1b_alignment);
        #print "@p1b_alignment\n"; print "$p1b_cov\n";
        $p1b_A_freq = 0;
        $p1b_C_freq = 0;
        $p1b_G_freq = 0;
        $p1b_T_freq = 0;
        if( $p1b_cov >= 5 ){
                for ( my $j=0; $j<=$#p1b_alignment; $j++ ) {
			if ( $p1b_alignment[$j] eq '-' or $p1b_alignment[$j] eq '+' ){
                                $j = $j + $p1b_alignment[$j+1] + 2;
                        } else {                        
			$p1b_A_freq++ if ($p1b_alignment[$j] =~ /[Aa]/ );
                        $p1b_C_freq++ if ($p1b_alignment[$j] =~ /[Cc]/ );
                        $p1b_G_freq++ if ($p1b_alignment[$j] =~ /[Gg]/ );
                        $p1b_T_freq++ if ($p1b_alignment[$j] =~ /[Tt]/ );
			}
                }
                #print "$p1b_A_freq\t$p1b_C_freq\t$p1b_G_freq\t$p1b_T_freq\n";
                $p1b_basecall = 'A' if ($p1b_A_freq/$p1b_cov > 0.9);
                $p1b_basecall = 'C' if ($p1b_C_freq/$p1b_cov > 0.9);
                $p1b_basecall = 'G' if ($p1b_G_freq/$p1b_cov > 0.9);
                $p1b_basecall = 'T' if ($p1b_T_freq/$p1b_cov > 0.9);
        }
        #print "$ref_name\t$ref_pos\t$p1b_basecall\n";
        $index1 = 1 if ( $p1a_basecall eq $p1b_basecall);
        $index1 = 2 if ( $p1a_basecall =~ /[ATCG]/ && $p1b_basecall eq 'X' );
        $index1 = 3 if ( $p1b_basecall =~ /[ATCG]/ && $p1a_basecall eq 'X' );
        $p1_mixbasecall = $p1a_basecall if ( $index1 == 1);
        $p1_mixbasecall = $p1a_basecall if ( $index1 == 2);
        $p1_mixbasecall = $p1b_basecall if ( $index1 == 3);
        #print "$p1a_basecall\t$p1b_basecall\t$p1_mixbasecall\n";
	if ( $p1_mixbasecall ne 'X' && $p1_mixbasecall ne $ref_base ){ 
        	print "$ref_name\t$ref_pos\t$ref_base\t$p1_mixbasecall\n";
	}


###########
=pod



### ir64_c_basecall ###
        $p2a_alignment =~ s/[\.\,]/$ref_base/g;
        @p2a_alignment = split('', $p2a_alignment);
        #print "@p2a_alignment\n"; print "$p2a_cov\n";
        $p2a_A_freq = 0;
        $p2a_C_freq = 0;
        $p2a_G_freq = 0;
        $p2a_T_freq = 0;
        if( $p2a_cov >= 10 ){
                for ( my $j=0; $j<=$#p2a_alignment; $j++ ) {
			if ( $p2a_alignment[$j] eq '-' or $p2a_alignment[$j] eq '+' ){
                                $j = $j + $p2a_alignment[$j+1] + 2;
                        } else {
                        $p2a_A_freq++ if ($p2a_alignment[$j] =~ /[Aa]/ );
                        $p2a_C_freq++ if ($p2a_alignment[$j] =~ /[Cc]/ );
                        $p2a_G_freq++ if ($p2a_alignment[$j] =~ /[Gg]/ );
                        $p2a_T_freq++ if ($p2a_alignment[$j] =~ /[Tt]/ );
			}
                }
                #print "$p2a_A_freq\t$p2a_C_freq\t$p2a_G_freq\t$p2a_T_freq\n";
                $p2a_basecall = 'A' if ($p2a_A_freq/$p2a_cov > 0.8);
                $p2a_basecall = 'C' if ($p2a_C_freq/$p2a_cov > 0.8);
                $p2a_basecall = 'G' if ($p2a_G_freq/$p2a_cov > 0.8);
                $p2a_basecall = 'T' if ($p2a_T_freq/$p2a_cov > 0.8);
        }



### ir64_s_basecall ###
        $p2b_alignment =~ s/[\.\,]/$ref_base/g;
        @p2b_alignment = split('', $p2b_alignment);
        #print "@p2b_alignment\n"; print "$p2b_cov\n";
        $p2b_A_freq = 0;
        $p2b_C_freq = 0;
        $p2b_G_freq = 0;
        $p2b_T_freq = 0;
        if( $p2b_cov >= 10 ){
                for ( my $j=0; $j<=$#p2b_alignment; $j++ ) {
			if ( $p2b_alignment[$j] eq '-' or $p2b_alignment[$j] eq '+' ){
                                $j = $j + $p2b_alignment[$j+1] + 2;
                        } else {
                        $p2b_A_freq++ if ($p2b_alignment[$j] =~ /[Aa]/ );
                        $p2b_C_freq++ if ($p2b_alignment[$j] =~ /[Cc]/ );
                        $p2b_G_freq++ if ($p2b_alignment[$j] =~ /[Gg]/ );
                        $p2b_T_freq++ if ($p2b_alignment[$j] =~ /[Tt]/ );
			}                
		}
                #print "$p1b_A_freq\t$p1b_C_freq\t$p1b_G_freq\t$p1b_T_freq\n";
                $p2b_basecall = 'A' if ($p2b_A_freq/$p2b_cov > 0.8);
                $p2b_basecall = 'C' if ($p2b_C_freq/$p2b_cov > 0.8);
                $p2b_basecall = 'G' if ($p2b_G_freq/$p2b_cov > 0.8);
                $p2b_basecall = 'T' if ($p2b_T_freq/$p2b_cov > 0.8);
        }
        $index2 = 1 if ( $p2a_basecall eq $p2b_basecall);
        $index2 = 2 if ( $p2a_basecall =~ /[ATCG]/ && $p2b_basecall eq 'X' );
        $index2 = 3 if ( $p2b_basecall =~ /[ATCG]/ && $p2a_basecall eq 'X' );
        $p2_mixbasecall = $p2a_basecall if ( $index2 == 1);
        $p2_mixbasecall = $p2a_basecall if ( $index2 == 2);
        $p2_mixbasecall = $p2b_basecall if ( $index2 == 3);
        #print "$p2a_basecall\t$p2b_basecall\t$p2_mixbasecall\n";
        
	
	
	
## basecall ##	
	if ( $p1_mixbasecall ne 'X' && $p2_mixbasecall ne 'X' ){
                if ( $p1_mixbasecall eq $p2_mixbasecall && $p1_mixbasecall ne $ref_base ){
			print "$ref_name\t$ref_pos\t$ref_base\t$p1_mixbasecall\t$p2_mixbasecall\n";
                }
        }
       

=cut

}
