#!/usr/bin/perl

use strict;
use warnings;

my $file = $ARGV[0];
my $ref_name;
my $ref_pos;
my $ref_base;
my $p1a_cov;
my $p1a_alignment;
my $p1b_cov;
my $p1b_alignment;
my $p2a_cov;
my $p2a_alignment;
my $p2b_cov;
my $p2b_alignment;


my $p1a_indel_freq;
my $p1b_indel_freq;
my $p2a_indel_freq;
my $p2b_indel_freq;

my @p1a_indel_array;
my @p1b_indel_array;
my @p2a_indel_array;
my @p2b_indel_array;




my %p1a_count;
my %p1b_count;
my %p2a_count;
my %p2b_count;

my $p1a_basecall;
my $p1b_basecall;
my $p2a_basecall;
my $p2b_basecall;

my @p1a_alignment;
my @p1b_alignment;
my @p2a_alignment;
my @p2b_alignment;

my $p1_mixbasecall;
my $p2_mixbasecall;

open IN, "<$file";



while (<IN>) {
        chomp;
	my $input = $_;
	$p1_mixbasecall = 'X';
	$p2_mixbasecall = 'X';

        my $index1 = 0;
        my $index2 = 0;
        if ($input =~ /(.+)\t(\d+)\t(\w)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)\t(\d+)\t(.+)\t(.+)/){
                $ref_name = $1;
                $ref_pos  = $2;
                $ref_base = $3;
                $p1a_cov  = $4;
                $p1a_alignment = $5;
                $p1b_cov = $7;
                $p1b_alignment = $8;
                $p2a_cov = $10;
                $p2a_alignment = $11;
                $p2b_cov = $13;
                $p2b_alignment = $14;
		
	}

### apo_c_basecall ###


        $p1a_alignment =~ s/[\.\,]/$ref_base/g;
        @p1a_alignment = split('', $p1a_alignment);

	$p1a_basecall = 'X';
	undef %p1a_count;
	@p1a_indel_array = ();
        $p1a_indel_freq = 0;

        if( $p1a_cov >= 5 ){
                for ( my $j=0; $j<=$#p1a_alignment; $j++  ) {
			if ( $p1a_alignment[$j] eq '-' or $p1a_alignment[$j] eq '+' ){
				my $indel_num = $p1a_alignment[$j+1];
				$indel_num = join '', @p1a_alignment[($j+1)..($j+2)] if( $p1a_alignment[$j+2] =~ /\d/ );
				my $indel = join '', @p1a_alignment[($j)..($j+2+$indel_num-1)];
				$indel = join '', @p1a_alignment[($j)..($j+3+$indel_num-1)] if( $p1a_alignment[$j+2] =~ /\d/ );	
				$indel = uc($indel);
				push(@p1a_indel_array, $indel);
				$j = $j + $p1a_alignment[$j+1] + 1;
			}elsif( $p1a_alignment[$j] eq '^' ){
	        		$j = $j + 1;
			}elsif( $p1a_alignment[$j] eq '$' ){
				$j = $j;
			}
		}
		if( $#p1a_indel_array >= 0 ){
			foreach my $seq (@p1a_indel_array) {
				$p1a_count{$seq}++;
			}

			my @p1a_indel_seq = reverse sort { $p1a_count{$a} <=> $p1a_count{$b} } keys %p1a_count; 
			$p1a_indel_freq = $p1a_count{$p1a_indel_seq[0]}/$p1a_cov;
                        $p1a_basecall = $p1a_indel_seq[0] if( $p1a_indel_freq > 0.5 );

		}

        }





### apo_s_basecall ###
        $p1b_alignment =~ s/[\.\,]/$ref_base/g;
        @p1b_alignment = split('', $p1b_alignment);

	undef %p1b_count;
	$p1b_basecall = 'X';
	@p1b_indel_array = ();
        $p1b_indel_freq = 0;

        if( $p1b_cov >= 5 ){
                for ( my $j=0; $j<=$#p1b_alignment; $j++  ) {
                        if ( $p1b_alignment[$j] eq '-' or $p1b_alignment[$j] eq '+' ){
                                my $indel_num = $p1b_alignment[$j+1];
                                $indel_num = join '', @p1b_alignment[($j+1)..($j+2)] if( $p1b_alignment[$j+2] =~ /\d/ );
                                my $indel = join '', @p1b_alignment[($j)..($j+2+$indel_num-1)];
                                $indel = join '', @p1b_alignment[($j)..($j+3+$indel_num-1)] if( $p1b_alignment[$j+2] =~ /\d/ );
                                $indel = uc($indel);
                                push(@p1b_indel_array, $indel);
                                $j = $j + $p1b_alignment[$j+1] + 1;
                        }elsif( $p1b_alignment[$j] eq '^' ){
                                $j = $j + 1;
                        }elsif( $p1b_alignment[$j] eq '$' ){
                                $j = $j;
                        }
                }
                if( $#p1b_indel_array >= 0 ){
                        foreach my $seq (@p1b_indel_array) {
                                $p1b_count{$seq}++;
                        }
                        my @p1b_indel_seq = reverse sort { $p1b_count{$a} <=> $p1b_count{$b} } keys %p1b_count;
			$p1b_indel_freq = $p1b_count{$p1b_indel_seq[0]}/$p1b_cov;
                        $p1b_basecall = $p1b_indel_seq[0] if( $p1b_indel_freq > 0.5 );
                }
        }



        $index1 = 1 if ( $p1a_basecall eq $p1b_basecall && $p1a_basecall ne 'X');
        $index1 = 2 if ( $p1a_basecall ne 'X' && $p1b_basecall eq 'X' );
	$index1 = 3 if ( $p1b_basecall ne 'X' && $p1a_basecall eq 'X' );
        $p1_mixbasecall = $p1a_basecall if ( $index1 == 1);
	$p1_mixbasecall = $p1a_basecall if ( $index1 == 2);
	$p1_mixbasecall = $p1b_basecall if ( $index1 == 3);




### ir64_c_basecall ###
        $p2a_alignment =~ s/[\.\,]/$ref_base/g;
        @p2a_alignment = split('', $p2a_alignment);

	undef %p2a_count;
	$p2a_basecall = 'X';
	@p2a_indel_array = ();
        $p2a_indel_freq = 0;

        if( $p2a_cov >= 5 ){
                for ( my $j=0; $j<=$#p2a_alignment; $j++  ) {
                        if ( $p2a_alignment[$j] eq '-' or $p2a_alignment[$j] eq '+' ){
                                my $indel_num = $p2a_alignment[$j+1];
                                $indel_num = join '', @p2a_alignment[($j+1)..($j+2)] if( $p2a_alignment[$j+2] =~ /\d/ );
                                my $indel = join '', @p2a_alignment[($j)..($j+2+$indel_num-1)];
                                $indel = join '', @p2a_alignment[($j)..($j+3+$indel_num-1)] if( $p2a_alignment[$j+2] =~ /\d/ );
                                $indel = uc($indel);
                                push(@p2a_indel_array, $indel);
                                $j = $j + $p2a_alignment[$j+1] + 1;
                        }elsif( $p2a_alignment[$j] eq '^' ){
                                $j = $j + 1;
                        }elsif( $p2a_alignment[$j] eq '$' ){
                                $j = $j;
                        }
                }
                if( $#p2a_indel_array >= 0 ){
                        foreach my $seq (@p2a_indel_array) {
                                $p2a_count{$seq}++;
                        }

                        my @p2a_indel_seq = reverse sort { $p2a_count{$a} <=> $p2a_count{$b} } keys %p2a_count;
			$p2a_indel_freq = $p2a_count{$p2a_indel_seq[0]}/$p2a_cov;
                        $p2a_basecall = $p2a_indel_seq[0] if( $p2a_indel_freq > 0.5 );
                }

        }



### ir64_s_basecall ###
        $p2b_alignment =~ s/[\.\,]/$ref_base/g;
        @p2b_alignment = split('', $p2b_alignment);

	undef %p2b_count;
	$p2b_basecall = 'X';
        @p2b_indel_array = ();
        $p2b_indel_freq = 0;

        if( $p2b_cov >= 5 ){
                for ( my $j=0; $j<=$#p2b_alignment; $j++  ) {
                        if ( $p2b_alignment[$j] eq '-' or $p2b_alignment[$j] eq '+' ){
                                my $indel_num = $p2b_alignment[$j+1];
                                $indel_num = join '', @p2b_alignment[($j+1)..($j+2)] if( $p2b_alignment[$j+2] =~ /\d/ );
                                my $indel = join '', @p2b_alignment[($j)..($j+2+$indel_num-1)];
                                $indel = join '', @p2b_alignment[($j)..($j+3+$indel_num-1)] if( $p2b_alignment[$j+2] =~ /\d/ );
                                $indel = uc($indel);
                                push(@p2b_indel_array, $indel);
                                $j = $j + $p2b_alignment[$j+1] + 1;
                        }elsif( $p2b_alignment[$j] eq '^' ){
                                $j = $j + 1;
                        }elsif( $p2b_alignment[$j] eq '$' ){
                                $j = $j;
                        }
                }
                if( $#p2b_indel_array >= 0 ){
                        foreach my $seq (@p2b_indel_array) {
                                $p2b_count{$seq}++;
                        }
                        my @p2b_indel_seq = reverse sort { $p2b_count{$a} <=> $p2b_count{$b} } keys %p2b_count;
                        $p2b_indel_freq = $p2b_count{$p2b_indel_seq[0]}/$p2b_cov;
			$p2b_basecall = $p2b_indel_seq[0] if( $p2b_indel_freq > 0.5 );
                        
                }
        }
	$index2 = 1 if ( $p2a_basecall eq $p2b_basecall && $p2a_basecall ne 'X');
        $index2 = 2 if ( $p2a_basecall ne 'X' && $p2b_basecall eq 'X' );
        $index2 = 3 if ( $p2b_basecall ne 'X' && $p2a_basecall eq 'X' );
        $p2_mixbasecall = $p2a_basecall if ( $index2 == 1);
        $p2_mixbasecall = $p2a_basecall if ( $index2 == 2);
        $p2_mixbasecall = $p2b_basecall if ( $index2 == 3);
        
	
	
	
## basecall ##	
	if ( $p1_mixbasecall ne 'X' && $p2_mixbasecall ne 'X' ){
                if ( $p1_mixbasecall eq $p2_mixbasecall ){
			print "$ref_name\t$ref_pos\t$p1_mixbasecall\n";
                }
        }
       



}
