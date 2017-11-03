use Bio::DB::Fasta;

&usage() unless scalar(@ARGV) == 2;

my $infile = $ARGV[0];
my $list = $ARGV[1];

my $db  = Bio::DB::Fasta->new( $infile );
#my $outseq = Bio::SeqIO->new(-file => ">$outfile", -format => "fasta");

open(INLIST, $list);

while(<INLIST>) {
  chomp;
  $SeqID = $_;
  my $FaSeq = $db->seq($SeqID);
    if  (!defined( $FaSeq )) {            die "Sequence $SeqID not found. \n"     }     print ">$SeqID\n", "$FaSeq\n";
}




sub usage {
print STDERR "
Extract fasta records using a list.

Usage:
  perl ExtractFastaRecords.pl fasta_file ID_list > output_file_name

where:
  fasta_file:  		the input fasta file
  ID_list:   		List of records to look for  
  output_file_name:  	name for the output file

#Note: Adapted from Stefano Berri's script found at https://www.biostars.org/p/2822/

";
exit(1);
}



