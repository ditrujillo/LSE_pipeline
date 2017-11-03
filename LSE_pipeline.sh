#!/bin/bash

# Author: Diana Trujillo ditrujillo@gmail.com

# MSI cluster researchers should load the following data:

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

# Note: recommended resources are 16 threads and 22 Gb mem

#########################################################################
################################# Setup ################################# 
#########################################################################

# Objective: Set up the data

# Prior to running Setup.sh: 
mkdir -p ~/LSE_pipeline/Setup/Mt
mkdir -p ~/LSE_pipeline/Setup/Gm
mkdir -p ~/LSE_pipeline/Setup/Tp
mkdir -p ~/LSE_pipeline/Setup/Pv
mkdir -p ~/LSE_pipeline/Setup/Ad

#Connect to JGI Phytozome to download genomes
#https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#
#Download to ~/LSE_pipeline/Setup/Mt: Mtruncatula_285_Mt4.0.fa.gz 
#Download to ~/LSE_pipeline/Setup/Gm: Gmax_275_v2.0.fa.gz
#Download to ~/LSE_pipeline/Setup/Tp: Tpratense_385_v2.fa.gz
#Download to ~/LSE_pipeline/Setup/Pv: Pvulgaris_442_v2.0.fa.gz

# Process:
# Run Setup.sh to do the following:
# Download RNAseq data for each species, and Arachis duranensis genome
# Build bowtie index for each genome
# Get FASTQC quality scores for the RNAseq and trim when necessary
# Map reads to the reference genome using Tophat
# Use Cufflinks to get nodule transcript boundaries

# Output:
# Nodule transcript boundaries: 
#   ~/LSE_pipeline/Setup/*/Nodules/Cufflinks/transcripts.gtf
# Mapped reads files from Ad,Pv,Gm,Mt,Tp (Species) X Flower,Leaf,Root,Nodule (Tissues):
#   ~/LSE_pipeline/Setup/*/*/accepted_hits.sorted.bam

cd ~/LSE_pipeline
./Setup.sh


#########################################################################
####################### Stage1: Nod-Specific SSPs ####################### 
#########################################################################

# Objective: Get a candidate pool of nodule-specific small, putatively secreted peptides

################################################
########### Stage1: Gene Specificity ########### 
################################################

# Objective: Get a list of nodule-specific transcripts

# Input:
# Nodule transcript boundaries: 
#   ~/LSE_pipeline/Setup/*/Nodules/Cufflinks/transcripts.gtf
# Mapped reads files from Ad,Pv,Gm,Mt,Tp (Species) X Flower,Leaf,Root,Nodule (Tissues):
#   ~/LSE_pipeline/Setup/*/*/accepted_hits.sorted.bam

# Process:
# For each tissue, find gene expression based on nodule transcriptome
# For each gene, get nodule 'relative expression' estimate --> nodule / (average of Flower, Leaf, Root)
# If the nodule 'relative expression' estimate is a ratio greater than 4 AND nodule FPKM is greater
#   than 10, the nodule transcript is considered 'nodule-specific'

# Output:
# List of nodule-specific transcripts for Ad, Pv, Gm, Mt, Tp:
#   ~/LSE_pipeline/Setup/*/NodSpecList.txt
# Fasta file with nodule-npecific transcripts:
#   ~/LSE_pipeline/Setup/*/*-NodSpec.fa


################################################
################ Stage1: Filter ################
################################################

# Objective: Within nodule-specific transcripts, find small, putatively secreted peptides 

# Input:
# List of nodule-specific transcripts for Ad,Pv,Gm,Mt,Tp:
#   ~/LSE_pipeline/Setup/*/NodSpecList.txt
# Fasta file with nodule-npecific transcripts:
#   ~/LSE_pipeline/Setup/*/*-NodSpec.fa

# Process:
# Translate DNA transcripts into all 6 frames
# Get all possible open reading frames (ORFs) from the nodule-specific transcripts
# Trim ORFs until the first methionine (M) residue
# Retain peptides between 35 to 250 amino acids in length
#   --> Longest known signaling peptides is systemin = 200 aa
#   --> 250 upper limit was chosen as a conservative size
# Find peptides that contain a signal peptide
# Some of the peptides within a single species have 100% identitiy. Retain unique peptides

# Output:
# Mature, small peptides from nodule-specific transcripts:
#   ~/LSE_pipeline/Stage1/SigP/*-mature

./Stage1.sh


#########################################################################
################ Stage2: Expand Pool of candidate genes #################
#########################################################################

# Objective: Find lineage-specific expansions (LSEs) of nodule-specific SSPs 
# and expand the pool of candidate genes within those LSEs.

# Process:
# Based on an initial pool of nodule-specific SSPs, cluster transcripts into 
# gene families based on sequence homology. Restrict gene families to those 
# with a LSE in a subset of species. Expand the list of candidate genes within
# a phylogeny-restricted family by doing homology searches against 
# reference genomes. Do this process iteratively to refine the homology search 
# outcome.

################################################
################ Stage2: Expand1 ###############
################################################

# Input:
# Mature, small peptides from nodule specific transcripts
#   ~/LSE_pipeline/Stage1/SigP/*-mature

# Process:
# Find homology between putative nodule-specific SSPs by conducting All-vs.-All BLAST
# Cluster genes based on sequence homology using MCL 
#   Inflation parameter set to 1.5 (MCL i=1.5) to cluster more loosely than default (i=2)
# Within MCL clusters of nodule-specific SSPs, count occurrences within each group for each 
#   species to get a 'Total Expansion Estimate'. When a cluster has 3 times or higher members 
#   in a species relative to another species, this cluster passes to the next step as a 
#   potential LSE
# Genes within potential LSE groups are aligned using ClustalW, then multiple sequence 
#   alignments are used to conduct homology-based searches against genomes using the Small  
#   Peptide Alignment Detection Application (SPADA)
# SPADA hits are filtered prior to going to the next stage, based on the following criteria:
#   - Nodule-specific
#   - 35 to 250 amino acids in length
#   - Contains signal peptide

# Output:
# Nodule-specific small secreted peptides from  gene families that have undergone Lineage-
#   Specific Expansions:
#   ~/LSE_pipeline/Stage2/Expand1/Spada/FilterSpadaHits/SpadaHitsNodSize.fa

### Note: This step relies on the Small Peptide Alignment Detection Application (SPADA) ###
#Prior to running SPADA, the SPADA configuration file ~/LSE_pipeline/ExampleFiles/Cluster.txt 
#must be set up to run on the user's computer, depending on the location of the program dependencies
#Stage2-Expand1.sh must be modified with the locations of the user files, in order to run spada.pl:

./Stage2-Expand1.sh


################################################
################ Stage2: Expand2 ###############
################################################

# Input:
# Nodule-specific small secreted peptides from  gene families that have undergone LSEs:
#   ~/LSE_pipeline/Stage2/Expand1/Spada/FilterSpadaHits/SpadaHitsNodSize.fa

# Process:
# Repeat the steps outlined in Stage2: Expand1, in order to expand the pool of genes for
#   candidate LSE families.
# The only difference with Stage2: Expand1, is a more restrictive inflation parameter during
#   MCL clustering (i=1.8)

# Output:
# Refined set of Nodule-specific small secreted peptides from  gene families that have 
# undergone LSEs:
#   ~/LSE_pipeline/Stage2/Expand2/Spada/FilterSpadaHits/SpadaHits.fa

#Note: This step relies on the Small Peptide Alignment Detection Application (SPADA)
#Stage2-Expand2.sh must be modified with the locations of the user files, in order to run spada.pl:

./Stage2-Expand2.sh


#########################################################################
############## Stage3: Curate candidate LSE gene families ###############
#########################################################################

# Manual annotation and curation of peptide clusters:
# In Stage2: Expand2, the peptides detected by SPADA were BLAST searched against M. truncatula 
#   gene annotations to assign a putative family ID to the MCL clusters:
#   ~/LSE_pipeline/Stage2/Expand2/Groups/Counts-TotExp-IDs.txt
# In addition to these annotations, clusters without a clear ID were searched against the
#   NCBI nr database using tblastn. In this way, peptides that actually corresponded to out of 
#   frame predictions of annotated genes were marked as 'DISCARD'
# The results of this manual annotation were compiled in the text file:
#   ~/LSE_pipeline/Stage2/Expand2/Groups/ExpansionAnalysis.txt
# An example file is provided: ExpansionAnalysis-Example.txt

# Manually Align Sequences:
# Gene families in ExpansionAnalysis.txt need to be manually aligned and visualized using Jalview 
# Alignments are done manually in order to retain the best quality sequences for a final HMM-based 
#   search of genomes using SPADA
# Curated gene family alignments are placed in:
#   ~/LSE_pipeline/Stage3/Seqs/Curated

# The commands in the following script need to be run on the command line interactively, in step 
# with the manual curation steps described above:

./Stage3.sh


# Finally, phylogenetic and expression data are combined to confirm LSEs:

./Heatmap.R




























































#Separate the expression data into groups
#After looking at expression profiles of detected genes, the following groups were dropped: CarbonicAnhydrase,Chitinase,Endogluc,HeatShockProt,Lectin,LectinLarge,LipidTransfer,LRR,LRRlarge,Pectinesterase,PhosphateResp,RipeningRel,XETHydrolase
cd ~/Data/Cluster/MCL2/Spada1/Expression/
mkdir -p Rfigs
for j in AeschNCRlike BowmanBirk CAPcysrich CysProtInhib DefensinRel ExtensinPRP GmNods LEED..PEED Leginsulin MBP1 NCRs NodGRP1 NodGRP2 NodPDPs nodGCRP PapainCystProt; do 
  awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}' > Rfigs/$j.tbl
  for i in Gm Mt Ad Tp Pv; do 
    grep "$j-" $i-summary.txt.ALL >> Rfigs/$j.tbl
    awk -F "\t" -v Sp=$i -v Grp=$j '{if (FNR==1) {next}; gsub("*", "", $19); h[$17]+=1 ; \
       if ($17 ~ Grp ) {print ">"Sp"|"$17"-"h[$17]"\n" $19}}' $i/61_final.gtb >> Rfigs/$j.fa
done; done

for j in AeschNCRlike BowmanBirk CAPcysrich CysProtInhib DefensinRel ExtensinPRP GmNods LEED..PEED Leginsulin MBP1 NCRs NodGRP1 NodGRP2 NodPDPs nodGCRP PapainCystProt; do 
  for i in Gm Mt Ad Tp Pv; do 
    awk -F "\t" -v Sp=$i -v Grp=$j '{if (FNR==1) {next}; gsub("*", "", $19); h[$17]+=1 ; \
       if ($17 ~ Grp ) {print ">"Sp"|"$17"-"h[$17]"\n" $19}}' $i/61_final.gtb >> Rfigs/$j.fa
done; done










#Trinity
#Species  Trinity   ORFs	35-350.uniq   SigP.uniq
#Tp	  290,445   5,421,911	692,826	      28,011
#Mt	  117,725   3,349,380	475,355       20,169
#Pv	   78,871     772,837	182,043        7,200
#Gm	   70,229     589,435	131,238        4,719
#Ad	  312,825   6,640,615	910,181       38,355

#Cufflinks and NodSpecList
#Species  Cufflnks  Nod-Spec	ORFs	  35-350.uniq   SigP.uniq
#Tp 	  62,559    4,655	133,733	  23,416	1,310
#Mt 	  39,344    2,360	127,523	  22,664	1,438
#Pv 	  45,579    1,837	 60,298	  13,557	  636
#Gm 	  48,694    3,488	 73,582	  17,125	  816
#Ad 	  40,316    1,738	 62,137	  13,452	  749

#Final Cluster Pipeline
#Species  Cufflnks  Nod-Spec	ORFs	  M35-250   SigP.uniq	Expand1
#Tp	  53,144    5,481	146,595   26,250    1,452	408
#Mt	  46,376    4,897	102,689   19,411    1,247	538
#Pv	  43,469    2,396	 74,918   16,761      775	 29
#Gm	  51,506    4,633	110,943   25,671    1,242	 64
#Ad	  40,316    2,114	 78,938   17,141      914	174

####Prev

#Expand1/Spada   #Prev		#New	#NewEval1	#NewEval10
#Ad/SpadaHits.fa:1119		1275	803
#Gm/SpadaHits.fa:1759		1909	1024
#Mt/SpadaHits.fa:1565		1798	1009
#Pv/SpadaHits.fa:912		1026	638
#Tp/SpadaHits.fa:1083		1238	814
#FPKM>0 AND 35-350		#New	#NewEval1	#NewEval10
#Ad-35-350.fa:763		789	651		666
#Gm-35-350.fa:935		982	756		781
#Mt-35-350.fa:983		1093	835		881
#Pv-35-350.fa:583		639	472		496
#Tp-35-350.fa:604		667	687		733

#awk 'BEGIN{while((getline<"groups.ProtIDs.txt")>0)l[">"$1]=1}/^>/{f=l[$1]?1:0}f' ../../goodProteins.fasta > groups.Prot.fa
#sort -k8,8nr -k1.6,1n counts-TotExp.txt | cut -f1 | head -80 > ExpansionList.txt
#module load t-coffee; cd ~/Data/Cluster/Expand1; mkdir -p tcoffee; cd tcoffee; 
#for i in $(cat ../Groups/ExpansionList.txt); do t_coffee ../Seqs/$i.fa -mode accurate; done


