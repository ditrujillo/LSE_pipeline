#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $fastafile = $ARGV[0];
my $count = $ARGV[1];
my $SplitInto = $ARGV[2];

#Obtain a list of the number of members for each group from counts.txt file
open(INCOUNT, $count) || die "Can't open fasta file '$fastafile'\n";

my @b = []; 
while(<INCOUNT>) {
    chomp;
    $_ =~ /(\d+)$/;
    push (@b, $1);
}

#Open the file I want to split and start counters
my $seqio = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");
my $groupNum = 1;
my $outseq = Bio::SeqIO->new(-file => ">Group$groupNum.fa", -format => "fasta");
my $groupMem = 1;

#Split the file
while(my $seqs = $seqio->next_seq) {

  unless ($groupNum > $SplitInto){
    if($groupMem > @b[$groupNum]) {
      $groupNum +=1;
      $groupMem = 1;
      $outseq = Bio::SeqIO->new(-file => ">Group$groupNum.fa", -format => "fasta") unless ($groupNum > $SplitInto);
      $outseq->write_seq($seqs) unless ($groupNum > $SplitInto);
      $groupMem +=1;
    }
    else {
      $outseq->write_seq($seqs);
      $groupMem +=1;
    }
  }
}

print @b[$groupNum];


sub usage {
print STDERR "
Based on counts.txt table, split fasta file into corresponding groups

Usage:
  perl SplitFastaByNum.pl fasta_file counts_table num

where:
  fasta_file:  		the input fasta file
  counts_table:		table in which the last column contains the number of fasta sequences per group
  num:			number of lines from counts_table to process

Output:
  fasta_files
  
";
exit(1);
}

