
#!/usr/bin/perl


use strict;
use warnings;
use Bio::SeqIO;

die "Usage: perl Temp.pl [Contigs FASTA]\n" unless @ARGV;

my ($ref_file, $basecall) = @ARGV;
my %ref_seq = ();


open IN2, "<$basecall";



### Read sequence in FASTA of both reference and contigs
my $seqio_obj = Bio::SeqIO -> new(-file => "$ref_file", -format => "fasta");
while (my $seq_obj = $seqio_obj -> next_seq) {
        my $id = $seq_obj -> display_id;
        my $seq = $seq_obj -> seq;
        die "The contigs contain irregular symbols in $id!\n" if ($seq =~ /\W/);
        $ref_seq{$id} = "$seq";
}



my @ref_name;
my @pos;
my @base_apo;


while( <IN2> ){
	chomp;
	my $line = $_;
	if ( $line =~ /(.+)\t(\d+)\t(\w)\t(\w)/){
		push(@ref_name, $1);
		push(@pos, $2);
		push(@base_apo, $4);
	}
}
close (IN2);

#substr($ref_seq{$ref_name[1]}, $pos[1]+1, 1) = $base_apo[1];
#substr($ref_seq{$ref_name[2]}, $pos[2]+1, 1) = $base_apo[2];
#substr($ref_seq{$ref_name[100092]}, $pos[100092]+1, 1) = $base_apo[100092];











my $i;
for ( $i=0; $i <=$#ref_name; $i++){
	substr($ref_seq{$ref_name[$i]}, $pos[$i]-1, 1) = $base_apo[$i];
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
