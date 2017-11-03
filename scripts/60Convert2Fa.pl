use strict;
use Bio::SeqIO;


my $file   = $ARGV[0];
my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.60", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
  my $id  = $seqs->display_id;
  my $seq = $seqs->seq;
    $outseq->write_seq($seqs);
}


sub usage {
print STDERR "
Convert fasta file to 60 characters per line

Usage:
  perl 60Convert2Fa.pl fasta_file

where:
  fasta_file:  		the input fasta file

Output:
  fasta_file.60
  
";
exit(1);
}
