#! /usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# Iker Irisarri, Jul 2019. MNCN-CSIC

# save HGs into hash

open(IN1, "<", "../2_PAPS_noisy/Output/Diapsida-atleast1_Mammalia-atleast1_outgroup-absent_3865_HGs_MCL_genes_IDs.out") or die "Can't open file!\n";

my %genes_HGs;

while ( my $line =<IN1> ) {
    chomp $line;
    my @lines = split ("\t", $line);
	
	my $HG = shift @lines; # get first element of line (HG)
	
	# save rest of elements
	foreach my $elem (@lines) {
		$genes_HGs{$elem} = $HG;
	}
}

#print Dumper \%genes_HGs;

# save queries (false positives) into hash

my %queries;

open(IN2, "<", "Amniota_novel_vs_nr_wo_used_genera.outfmt6.allqueries") or die "Can't open file!\n";
while ( my $q =<IN2> ) {
    chomp $q;
	$queries{$q} = "1";
}


# print out HGs where queries are found

my %output;

foreach my $query ( keys %queries ) {

	if ( exists $genes_HGs{$query} ) {
	
		my $hg = $genes_HGs{$query};
		
		print "$hg\n"; # print output
		
		$output{$hg}{$query} = 1; # save for extended output
	}
}

# generate outfile with extended output

open(OUT, ">", "Amniota_novel_HG_with_false_positives.txt");

foreach my $k ( sort keys %output ) {

	print OUT "$k";
	
	my %inner = %{ $output{$k} }; # derreference

	foreach my $m ( sort keys %inner ) {
		
		print OUT "\t$m";
	}
	print OUT "\n";
}