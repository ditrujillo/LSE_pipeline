#!/usr/bin/perl

use strict;    
use warnings;  


open (my $handle, "<@ARGV");

my $header = "";
my $sequence = "";
my $seq2 = "";
my $count = 0;
my $line = "";

while (<$handle>) {

    $line = $_;
    chomp($line);

    if(substr($line,0,1) eq '>') {
	$header = $line;
	$seq2 = $sequence;
	$seq2 =~ s/^[A-L,N-Z]*M/M/;
	for (my $i=0; $i<length($seq2); $i+=60) {
		print substr($seq2,$i,60);
		print "\n";
	}
	print $header."\n";
	$sequence = "";	
    } else {	
	$sequence .= $line;
    }
}

$seq2 = $sequence;
$seq2 =~ s/^[A-L,N-Z]*M/M/;
for (my $i=0; $i<length($seq2); $i+=60) {
	print substr($seq2,$i,60);
	print "\n";
}


close ($handle);

# 	$count ++;
#print "Number of entries: $count \n";


sub usage {
print STDERR "
Trim amino acid sequence until the first M residue is detected. Do nothing if M is not detected at all

Usage:
  perl TrimUntilM.pl ORF_fasta > Trimmed_ORF_fasta

where:
  ORF_fasta:  		fasta file name with all possible ORFs
  Trimmed_ORF_fasta:	output file with trimmed ORFs


";
exit(1);
}



