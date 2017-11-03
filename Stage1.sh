#!/bin/bash -l

module load fastqc/0.11.5
module load sratoolkit/2.8.2
module load bowtie/2.3.0
module load tophat/2.0.13
module load samtools/1.3
module load cufflinks/2.0.0
module load ncbi_blast+/2.2.28
module load mcl/10.201 
module load bioperl/1.6.901
module load kent/2.4.7
module load clustalw/2.1


#########################################################################
####################### Stage1: Nod-Specific SSPs ####################### 
#########################################################################

################################################
########### Stage1: Gene Specificity ########### 
################################################

#Get list of transcripts with nodule-specific expression

for j in Ad Mt Tp Gm Pv; do

  #For each tissue, find gene expression based on nodule transcriptome
  cd ~/LSE_pipeline/Setup/$j
  for i in Flower Leaf Root Nodules; do 
  cufflinks -u --min-intron-length 20 --max-intron-length 2000 -G Nodules/Cufflinks/transcripts.gtf -p 16 -o $i/CufflinksNodGuided \
      --max-bundle-frags 50000000 $i/accepted_hits.sorted.bam
  done

  #Extract nodule transcripts
  gffread -g $j.fa -w Nodules/Cufflinks/$j\Exons.fa Nodules/Cufflinks/transcripts.gtf

  #For each gene, get nodule 'relative expression' estimate --> nodule / (average of Flower, Leaf, Root)
  tail -n +2 Nodules/CufflinksNodGuided/genes.fpkm_tracking  | cut -f1,7 > summary.txt
  for i in Flower Leaf Root Nodules; do 
    awk -v OFS="\t" 'NR==FNR {h[$7] = $10; next} {print $0,h[$2]}' $i/CufflinksNodGuided/genes.fpkm_tracking summary.txt > summary.$i
    mv summary.$i summary.txt
  done

  #Get list of 'nodule-specific' transcripts for Ad (2122), Pv (2426), Gm (4305), Mt (4951), Tp (5602):
  awk -v OFS="\t" 'BEGIN{print "Protein\tLocus\tFlower\tLeaf\tRoot\tNodule\tAvgNon\tRelExp\tSpecificity"}
     {Spec=0; avgnon=($3+$4+$5)/3; relexp=$6/(avgnon+0.1)
     if ($6 > 10 && relexp > 4) {Spec=1}
     print $0, avgnon, relexp, Spec 
  }' summary.txt > TotExp.txt
  awk '$9==1{print $1".1"}' TotExp.txt > NodSpecList.txt
  perl ../../scripts/ExtractFastaRecords.pl Nodules/Cufflinks/$j\Exons.fa NodSpecList.txt > $j-NodSpec.fa

done


################################################
################ Stage1: Filter ################
################################################

#Symlink copy Nod-specific transcripts
mkdir -p ~/LSE_pipeline/Stage1/NodTranscripts
cd ~/LSE_pipeline/Stage1/NodTranscripts
cp -s ~/LSE_pipeline/Setup/*/*-NodSpec.fa .

#Translate DNA transcripts: translate6.pl
#Get all possible open reading frames (ORFs): getORFs.pl
mkdir -p ../ORFs; cd ../ORFs
for i in Tp Mt Pv Gm Ad; do
  perl ../../scripts/translate6.pl ../NodTranscripts/$i-NodSpec.fa $i-6.fa 
  perl ../../scripts/getORFs.pl $i-6.fa $i-ORFs.fa
done

#Trim ORFs until first M residue: TrimUntilM.pl
#Retain peptides with length of 35-250: faFilter
#Retain unique peptides: Collapse100.pl
mkdir -p ../TrimMSize; cd ../TrimMSize
for i in Tp Mt Pv Gm Ad; do
  perl ../../scripts/TrimUntilM.pl ../ORFs/$i-ORFs.fa > $i-ORFsM.fa
  faFilter -minSize=35 -maxSize=250 $i-ORFsM.fa $i-ORFsMSize.fa 		 
  perl ../../scripts/Collapse100.pl $i-ORFsMSize.fa
  #Merge into one file
  awk -v Sp=$i '{gsub(/\|/,"_", $0); gsub(/>/,">"Sp"|", $0); {print $0}}' $i-ORFsMSize.fa.uniq >> ORFsMSize.fa
done 

#Find peptides that contain a signal peptide: splitJobSigP.pl
#Retain unique peptides: Collapse100.pl
mkdir -p ../SigP
cd ../SigP
for i in Tp Mt Pv Gm Ad; do
  mkdir -p $i; cd $i
  perl ../../../scripts/splitJobSigP.pl -p signalp4 -n 1000 ../../TrimMSize/$i-ORFsMSize.fa.uniq 
  cat mature* > ../$i-mature
  cd ..; rm -rf $i
  perl ../../../scripts/Collapse100.pl $i-mature
done



