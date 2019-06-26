#!/usr/bin/perl -w
use strict;
use warnings;
###
#set file input from command line argument
my $sample = $ARGV[0];
my $dir = $ARGV[1];
my $genesfile = $ARGV[2];
my $file = "$dir/splicing/ReadsPerGene.out.tab";
###
my $strandedoutput = "$dir/genecounts/$sample.nodups.274genes.ncbi.txt";
###
open my $countsfile, '<', $file or die "Can't read $file : $!";
my @counts = <$countsfile>;
close $countsfile;
open my $geneslist, '<', $genesfile or die "Can't read $genesfile : $!";
my @genes = <$geneslist>;
close $geneslist;
###
open my $strandedfile, '>', $strandedoutput or die "Can't read $strandedoutput : $!";
###
COUNTS: foreach my $countline (@counts){
	chomp $countline;
	my ($gene, $unstranded, $na, $stranded) = split('\t', $countline);
	GENE: foreach my $currentgene (@genes){
		chomp $currentgene;
		my ($geneid, $ensgid, $genelength) = split('\t', $currentgene);
		if ($gene eq $geneid){
			my $strandedline = join("\t", $geneid, $stranded);
			print $strandedfile "$strandedline\n";
			next COUNTS;
		} else {
			next GENE;
		}
	}
}
close $strandedfile;
