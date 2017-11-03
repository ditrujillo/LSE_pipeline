use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;

my $fi   = $ARGV[0];
my $fo   = $ARGV[1];

my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');


my $cutoff_missing ||= 0.5;
my $sep ||= "|";

while(my $seq = $seqHI->next_seq()) {

  my ($id, $seqStr) = ($seq->id, $seq->seq);
  my ($beg, $end, $srd) = (1, 3*length($seqStr), "+");

        if($seq->id =~ /^(\S+)\Q$sep\E(\d+)\Q$sep\E(\d+)\Q$sep\E([\+\-])$/) {
            ($id, $beg, $end, $srd) = ($1, $2, $3, $4);
        }
        while( $seqStr =~ /([^\*]{15,})/g ) {
            my ($begL, $endL) = ($-[1]+1, $+[1]);
            my ($begG, $endG);
            if($srd eq "-") {
                $begG = $end - $endL*3 + 1;
                $endG = $end - ($begL*3-2) + 1;
            } else {
                $begG = $beg + ($begL*3-2) - 1;
                $endG = $beg + $endL*3 - 1;
            }
            my $sseq = $1;
            my $sid = join($sep, $id, "$begG-$endG", $srd, "x");
            my $n_x =()= $1 =~ /X/gi;
            if($n_x / ($endL-$begL+1) <= $cutoff_missing) {
                $seqHO->write_seq(Bio::Seq->new(-id=>$sid, -seq=>$sseq));
            }
        }
    }
    $seqHI->close();
    $seqHO->close();


sub usage {
print STDERR "
Get all possible open reading frames (ORFs) from a translated DNA sequence

Usage:
  perl getORFs.pl 6_frame_fasta ORF_fasta

where:
  6_frame_fasta:  	the input amino acid fasta file of 6 frame translated DNA
  ORF_fasta:  		output file name with all possible ORFs

#Note: this script is adapted from Peng Zhou's PrepareGenome.pm package within SPADA 
#Reference: Zhou, Peng et al. “Detecting Small Plant Peptides Using SPADA (Small Peptide
#              Alignment Discovery Application).” BMC bioinformatics 14 (2013): 335 

";
exit(1);
}


