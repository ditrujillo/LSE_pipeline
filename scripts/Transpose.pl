#!/usr/bin/perl

use strict;    
use warnings;  

open(my $GROUPS, "<$ARGV[0]");
open(my $FUNCTION, "<$ARGV[1]");

my @list1=<$GROUPS>;
my @list2=<$FUNCTION>;

close($GROUPS);
close($FUNCTION);

chomp(@list1,@list2);

my $linecount = 0;

foreach my $line1 (@list1) {

    my $NA = "NA";
    my @array = split(/\s+/, $line1);
    $line1 =~ /(\d+)$/;
    my $MemNum = $1;
    print $line1."\t"; #.$MemNum

    for (my $i=0; $i < $MemNum; $i++){
#    foreach my $line2 (@list2) {
#	my @array2 = split(/\s+|\:/, $line2);
#        if (grep /$array[1]/i,$line2){$groupID = $array2[0];}
#    }
#    print $array[0]."\t".$array[1]."\t".$groupID."\n";
	my @array2 = split(/\s+/, @list2[$linecount]);
	$linecount ++;
	if ($array2[1] eq ""){
           print "\t".$NA;
        }
        else {
           print "\t".$array2[1];
        }
    }
    print "\n";
}


sub usage {
print STDERR "
Transpose gene annotations to the counts.txt table

Usage:
  perl Transpose.pl counts_table function_table > counts_table_plus_functions

where:
  counts_table:  	table with number of genes per MCL cluster per species
  function_table:	for each transcript, putative annotation based on best BLAST hit


";
exit(1);
}



