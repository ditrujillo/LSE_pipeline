use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Log::Log4perl;
use Data::Dumper;
use List::Util qw/min max sum/;


my (%genetic_code) = (
  'TCA' => 'S',    # Serine
  'TCC' => 'S',    # Serine
  'TCG' => 'S',    # Serine
  'TCT' => 'S',    # Serine
  'TTC' => 'F',    # Phenylalanine
  'TTT' => 'F',    # Phenylalanine
  'TTA' => 'L',    # Leucine
  'TTG' => 'L',    # Leucine
  'TAC' => 'Y',    # Tyrosine
  'TAT' => 'Y',    # Tyrosine
  'TAA' => '*',    # Stop
  'TAG' => '*',    # Stop
  'TGC' => 'C',    # Cysteine
  'TGT' => 'C',    # Cysteine
  'TGA' => '*',    # Stop
  'TGG' => 'W',    # Tryptophan
  'CTA' => 'L',    # Leucine
  'CTC' => 'L',    # Leucine
  'CTG' => 'L',    # Leucine
  'CTT' => 'L',    # Leucine
  'CCA' => 'P',    # Proline
  'CCC' => 'P',    # Proline
  'CCG' => 'P',    # Proline
  'CCT' => 'P',    # Proline
  'CAC' => 'H',    # Histidine
  'CAT' => 'H',    # Histidine
  'CAA' => 'Q',    # Glutamine
  'CAG' => 'Q',    # Glutamine
  'CGA' => 'R',    # Arginine
  'CGC' => 'R',    # Arginine
  'CGG' => 'R',    # Arginine
  'CGT' => 'R',    # Arginine
  'ATA' => 'I',    # Isoleucine
  'ATC' => 'I',    # Isoleucine
  'ATT' => 'I',    # Isoleucine
  'ATG' => 'M',    # Methionine
  'ACA' => 'T',    # Threonine
  'ACC' => 'T',    # Threonine
  'ACG' => 'T',    # Threonine
  'ACT' => 'T',    # Threonine
  'AAC' => 'N',    # Asparagine
  'AAT' => 'N',    # Asparagine
  'AAA' => 'K',    # Lysine
  'AAG' => 'K',    # Lysine
  'AGC' => 'S',    # Serine
  'AGT' => 'S',    # Serine
  'AGA' => 'R',    # Arginine
  'AGG' => 'R',    # Arginine
  'GTA' => 'V',    # Valine
  'GTC' => 'V',    # Valine
  'GTG' => 'V',    # Valine
  'GTT' => 'V',    # Valine
  'GCA' => 'A',    # Alanine
  'GCC' => 'A',    # Alanine
  'GCG' => 'A',    # Alanine
  'GCT' => 'A',    # Alanine
  'GAC' => 'D',    # Aspartic Acid
  'GAT' => 'D',    # Aspartic Acid
  'GAA' => 'E',    # Glutamic Acid
  'GAG' => 'E',    # Glutamic Acid
  'GGA' => 'G',    # Glycine
  'GGC' => 'G',    # Glycine
  'GGG' => 'G',    # Glycine
  'GGT' => 'G',    # Glycine
  'NAA' => 'X',    
  'NAT' => 'X',    
  'NAC' => 'X',    
  'NAG' => 'X',    
  'NTA' => 'X',    
  'NTT' => 'X',    
  'NTC' => 'X',    
  'NTG' => 'X',    
  'NGA' => 'X',    
  'NGT' => 'X',    
  'NGC' => 'X',    
  'NGG' => 'X',    
  'NCA' => 'X',    
  'NCT' => 'X',    
  'NCC' => 'X',    
  'NCG' => 'X',    
  'ANA' => 'X',    
  'ANT' => 'X',    
  'ANG' => 'X',    
  'ANC' => 'X',    
  'TNA' => 'X',    
  'TNT' => 'X',    
  'TNG' => 'X',    
  'TNC' => 'X',    
  'GNA' => 'X',    
  'GNT' => 'X',    
  'GNG' => 'X',   
  'GNC' => 'X',    
  'CNA' => 'X',  
  'CNT' => 'X',   
  'CNG' => 'X',    
  'CNC' => 'X',    
  'AAN' => 'X',
  'ATN' => 'X', # in 3 of 4 codes for I otherwise M-start
  'AGN' => 'X',
  'ACN' => 'T',
  'TAN' => 'X',
  'TTN' => 'X',
  'TGN' => 'X',
  'TCN' => 'S',
  'GAN' => 'X',
  'GTN' => 'V',
  'GGN' => 'G',
  'GCN' => 'A',
  'CAN' => 'X',
  'CTN' => 'L',
  'CGN' => 'R',
  'CCN' => 'P',
  'ANN' => 'X',	
  'TNN' => 'X',	
  'GNN' => 'X',
  'CNN' => 'X',
  'NAN' => 'X',
  'NTN' => 'X',
  'NGN' => 'X',
  'NCN' => 'X',
  'NNA' => 'X',
  'NNT' => 'X',
  'NNG' => 'X',
  'NNC' => 'X',
  'NNN' => 'X',
);



my $fi   = $ARGV[0];
my $fo   = $ARGV[1];

my $seqHI = Bio::SeqIO->new(-file=>"<$fi", -format=>'fasta');
my $seqHO = Bio::SeqIO->new(-file=>">$fo", -format=>'fasta');

my $sep ||= "|";

while(my $seqO = $seqHI->next_seq()) {
    my ($id, $seq) = ($seqO->id, $seqO->seq);
    $id =~ s/\|/\_/;
    my $len = length($seq);

    for my $i (0..2) {
      my $n_codon = int( ($len - $i) / 3 );
      my $lenC = $n_codon * 3;
      my ($beg, $end) = ($i+1, $i+$lenC);
      my $dna = substr($seq, $beg-1, $lenC);
      my $pro = translate($dna);
      my $sid = join($sep, $id, $beg, $end, "+");
      $seqHO->write_seq( Bio::Seq->new(-id=>$sid, -seq=>$pro) );
    }
    $seq = revcom($seq);
    for my $i (0..2) {
      my $n_codon = int( ($len - $i) / 3 );
      my $lenC = $n_codon * 3;
      my ($beg, $end) = ($i+1, $i+$lenC);
      my $dna = substr($seq, $beg-1, $lenC);
      my $pro = translate($dna);
      my ($begO, $endO) = ($len-$end+1, $len-$beg+1);
      my $sid = join($sep, $id, $begO, $endO, "-");
      $seqHO->write_seq( Bio::Seq->new(-id=>$sid, -seq=>$pro) );
    }
}
$seqHI->close();
$seqHO->close();


sub revcom {
  my($dna) = @_;
  # Reverse the sequence
  my $revcom = reverse($dna);
  # Complement the sequence
  $revcom =~ tr/ACGTacgt/TGCAtgca/;
  return $revcom;
}

sub translate {
  my ($dna) = @_;
  my $seqlen = length($dna);
  my $n_codon = int($seqlen / 3);
  my $pro = '';
  for(my $i=0; $i < $n_codon ; $i ++) {
    my $codon = uc( substr($dna, $i*3, 3) );
    my $aa = exists $genetic_code{$codon} ? $genetic_code{$codon} : "X";
    $pro .= $aa;
  }
  return $pro;
}



sub usage {
print STDERR "
Convert fasta file to 60 characters per line

Usage:
  perl 60Convert2Fa.pl fasta_file
  perl translate6.pl fasta_file 6_frame_fasta

where:
  fasta_file:  		the input DNA fasta file
  6_frame_fasta:	the output amino acid fasta file of translated DNA to all 6 frames

#Note: this script is adapted from Peng Zhou's Seq.pm package within SPADA 
#Reference: Zhou, Peng et al. “Detecting Small Plant Peptides Using SPADA (Small Peptide
#              Alignment Discovery Application).” BMC bioinformatics 14 (2013): 335 

 
";
exit(1);
}





