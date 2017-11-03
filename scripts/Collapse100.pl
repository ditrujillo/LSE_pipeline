use strict;
use Bio::SeqIO;

my %unique;
my $file   = "@ARGV";
my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");
my $outseq2 = Bio::SeqIO->new(-file => ">$file.rep", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
  my $id  = $seqs->display_id;
  my $seq = $seqs->seq;
  if(exists($unique{$seq})) {
    $outseq2->write_seq($seqs);
  }
  else {
    $outseq->write_seq($seqs);
    $unique{$seq} +=1;
  }
}


sub usage {
print STDERR "

In a fasta file with multiple sequences, discard the sequences that have 100% identity with another

Usage:

  perl Collapse100.pl fasta_file

where:
  fasta_file:  	file with multiple fasta sequences

Output:
  fa.uniq:	fasta file with unique sequences
  fa.rep:	fasta file with repeated sequences

";
exit(1);
}


