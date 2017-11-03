#!/usr/bin/perl

use strict; 
use warnings;

# Script splits a job that takes as input a FASTA file
# into multiple smaller jobs with fewer input seqs/job.
#
# USAGE
#    splitJobIntoFewerInputSeqs.pl \
#	-p [targetp | signalp3 | signalp4 | phobius] \
#	-n <number of input seqs per job [Default 100]> \
#	<input FASTA file>
#
# Default values: 
#	-n 100
# --------------------------------------------------------------
use Getopt::Std;
use Bio::SeqIO;

# Directory where the split FASTA files ($numInputSeq per file)
# are temporarily stored.
my $tmpDir = ".splitJobIntoFewerInputSeqs";


# Get options
my(%opts,$key,$val);
getopts('n:p:', \%opts);
my $program;
my $numSeqsPerJob=100;
while(($key,$val) = each(%opts)) {
  if ($key eq "p") {$program = $val; next;}
  if ($key eq "n") {$numSeqsPerJob = $val; next;}
  print STDERR "Unrecognized option -$key"; die;
}
if (! defined $program) {die "Need to set option -p [targetp|signalp]:$!";}

my $command;
if ($program eq 'targetp') {
  $command = 'targetp -P';
} elsif ($program eq 'signalp3') {
  $command = 'signalp -t euk -f short -trunc 70';
} elsif ($program eq 'signalp4') {
  $command = 'signalp -t euk -m mature';

} elsif ($program eq 'phobius') {
  $command = 'phobius.pl -short';
} else {
  die "Need to set option -p [targetp|signalp3|signalp4]:$!";
}
print "Using command $command\n";

my $inputFASTAfile = shift @ARGV;
if (! -r $inputFASTAfile) {die "Need readable input FASTA file:$!";}
my $instream  = Bio::SeqIO->new(-file => "$inputFASTAfile", -format => 'Fasta');

# Create temporary FASTA files with numSeqsPerJob seqs per file
mkdir $tmpDir || die "Couldn't create temporary directory $tmpDir:$!";

my $fileNum = 0;
my $numReadSeqs = 0;
my $outstream;
while (my $seq = $instream->next_seq()) {
  # Time to start a new file?
  if ($numReadSeqs % $numSeqsPerJob == 0) {
    $fileNum++;
    my $outFile = $tmpDir . '/in.' . $fileNum . ".fa";
    $outstream = Bio::SeqIO->new(-file => ">$outFile", -format => 'Fasta');
  }

  $outstream->write_seq($seq);
  $numReadSeqs++;

}

# Now run $program
foreach my $infile (<$tmpDir/in.*.fa>) {
  if ($infile !~ /^$tmpDir\/in\.(\d+)\.fa$/) {die "How did this happen? $!";}
  my $fileNum = $1;

  my $outFile = "$program.$fileNum.out";
  my $comm = $command . "$outFile $infile > $outFile";  

  # Note that /bin/sh will obliterate $outFile if it exists.
  print "$comm\n";
  system ($comm);
}

# Cleanup: rm files and temp directory
#foreach my $infile (<$tmpDir/*>) { unlink $infile; }
#rmdir $tmpDir;  #BUG: this directory does not get removed, but it is empty


sub usage {
print STDERR "

Split a fasta file into multiple sequences, to run SignalP on s smaller number of sequences at a time

Usage:

  perl splitJobSigP.pl -p program -n num_seqs fasta_file

where:
  program:  	program that will be run
  num_seqs:	number of sequences in each temporary file
  fasta_file:	fasta file with amino acid sequences

Output:
  fa-mature:	fasta file with mature peptides

#Note: this script was written by Karen Tang - University of Minnesota 

";
exit(1);
}

