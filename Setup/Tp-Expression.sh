#!/bin/bash -l

# Author: Diana Trujillo ditrujillo@gmail.com

module load bowtie/2.3.0
module load sratoolkit/2.8.2
module load tophat/2.0.13
module load samtools/1.3
module load cufflinks/2.0.0
module load fastqc/0.11.5
module load bioperl/1.6.901


#Download Medicago truncatula 4.0 genome
mkdir -p ~/LSE_pipeline/Setup/Tp/fastqc
cd ~/LSE_pipeline/Setup/Tp
#Connect to JGI Phytozome 
#https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#
#Download to ~/LSE_pipeline/Setup/Tp: Tpratense_385_v2.fa.gz
gzip -d Tpratense_385_v2.fa.gz
mv Tpratense_385_v2.fa Tp.fa


#Build Bowtie Index
bowtie2-build Tp.fa Tp


#Download data
#Leaf 
for i in SRR2899578 SRR2899580 SRR2899625 SRR2899654 SRR2899603 SRR2899588 SRR2899632; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR289/$i/$i.sra
 fastq-dump $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
 cat $i.fastq >> Leaf.fastq
done

#Flower
for i in SRR2899677 SRR2899634 SRR2899626 SRR2899627; do 
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR289/$i/$i.sra
 fastq-dump $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
 cat $i.fastq >> Flower.fastq
done

#Root 
for i in SRR2899631 SRR2899633 SRR2899630 SRR2899628 SRR2899653; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR289/$i/$i.sra
 fastq-dump $i.sra
 fastqc $i.fastq; mv *fastqc.* fastqc/.
 cat $i.fastq >> Root.fastq
done

#Nodule
for i in SRR6251246; do
 wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR625/$i/$i.sra
 fastq-dump $i.sra --split-3; rm $i.sra
 fastqc $i\_1.fastq; mv *fastqc.* fastqc/.
 fastqc $i\_2.fastq; mv *fastqc.* fastqc/.
 mv $i\_1.fastq Nodule_1.fastq
 mv $i\_2.fastq Nodule_2.fastq
done


#Fastqc reports suggest a lot of adapter contamination for Flower, Leaf, Root data
#Trim:
java -jar ~/software/Trimmomatic/trimmomatic.jar SE -phred33 Leaf.fastq  Leaf_Trim.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70 ILLUMINACLIP:/home/youngn/truji033/software/Trimmomatic/adapters/TruSeq-ALL.fa:2:30:10
fastqc Leaf_Trim.fastq
java -jar ~/software/Trimmomatic/trimmomatic.jar SE -phred33 Flower.fastq  Flower_Trim.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70 ILLUMINACLIP:/home/youngn/truji033/software/Trimmomatic/adapters/TruSeq-ALL.fa:2:30:10
fastqc Flower_Trim.fastq
java -jar ~/software/Trimmomatic/trimmomatic.jar SE -phred33 Root.fastq  Root_Trim.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70 ILLUMINACLIP:/home/youngn/truji033/software/Trimmomatic/adapters/TruSeq-ALL.fa:2:30:10
fastqc Root_Trim.fastq
mv *fastqc.* fastqc/.

rm SRR213* SRR214* SRR625* 

#Tophat
tophat -o Leaf -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 24 --segment-mismatches 1 Tp Leaf_Trim.fastq
tophat -o Flower -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 24 --segment-mismatches 1 Tp Flower_Trim.fastq
tophat -o Root -i 20 -I 2000 --b2-very-sensitive -p 16 --segment-length 24 --segment-mismatches 1 Tp Root_Trim.fastq
tophat -o Nodules -i 20 -I 2000 --b2-very-sensitive -p 16 Tp Nodule_1.fastq Nodule_2.fastq


#Cufflinks
for i in Nodules Root Flower Leaf; do
 cd $i
 samtools sort accepted_hits.bam -o accepted_hits.sorted.bam
 samtools index accepted_hits.sorted.bam
 rm accepted_hits.bam unmapped.bam
 cufflinks -u --min-intron-length 20 --max-intron-length 2000 -p 16 -o Cufflinks --max-bundle-frags 50000000 \
   accepted_hits.sorted.bam
 cd ..
done


