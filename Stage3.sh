#!/bin/bash -l

# Author: Diana Trujillo ditrujillo@gmail.com

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
############## Stage3: Curate candidate LSE gene families ###############
#########################################################################

################ Manual Cluster ################

# Manual annotation and curation of peptide clusters:
# In Stage2: Expand2, the peptides detected by SPADA were BLAST searched against M. truncatula 
#   gene annotations to assign a putative family ID to the MCL clusters:
#   ~/LSE_pipeline/Stage2/Expand2/Groups/Counts-TotExp-IDs.txt
# In addition to these annotations, clusters without a clear ID were searched against the
#   NCBI nr database using tblastn. In this way, peptides that actually corresponded to out of 
#   frame predictions of annotated genes were marked as 'DISCARD'
# The results of this manual annotation were compiled in the text file:
#   ~/LSE_pipeline/Stage2/Expand2/Groups/ExpansionAnalysis.txt

mkdir -p ~/LSE_pipeline/Stage3/Seqs/Raw ~/LSE_pipeline/Stage3/Seqs/Curated ~/LSE_pipeline/Stage3/Seqs/Discarded
cd ~/LSE_pipeline/Stage3/Seqs
cp ~/LSE_pipeline/Stage2/Expand2/Spada/FilterSpadaHits/SpadaHits.fa .
cp ~/LSE_pipeline/Stage2/Expand2/Groups/ExpansionAnalysis.txt .

#Spada Hits grouped based on annotation in ExpansionAnalysis.txt
#An example file is provided: ExpansionAnalysis-Example.txt
while IFS='|' read Grp Family ; do 
 for i in $Grp; do
   ~/DNC/bioawk/bioawk -c fastx -v grp=$i '{if ($name~grp"-") {print ">"$name; print $seq}}' SpadaHits.fa >> Raw/$Family.fa; done
done < ExpansionAnalysis.txt

#24 different groups + 'DISCARD':
#AeschNCRlike, BowBirk, CAMlike, CAPCysRich, ChiHevWheatin, Cystatin, Defensin, DISCARD, EarlyNodulin, emp24, GlycoHydro, HeatShockProt, Kunitz, LEEDPEED, Leginsulin, LipidTransferProt(LTP), MBP1, NCRs, NodGRP, PDP, Peroxidase, PhaseoleaeNodulin, PRP, RipeningRelated, Thioredoxin


################# Manual Align #################

#Manually Align Sequences:
# 24 gene families in ~/LSE_pipeline/Stage3/Seqs/Raw were manually aligned and visualized using Jalview 
# Alignments were done manually in order to retain the best quality sequences for a final HMM-based 
#   search of genomes using SPADA
# Curated gene family alignments were placed in:
#   ~/LSE_pipeline/Stage3/Seqs/Curated
# Gene families in which most members were too long (>250 amino acids), or in which the nodule-specific 
#   members correspoded to truncated versions of larger proteins, were placed in 
#   ~/LSE_pipeline/Stage3/Seqs/Discarded

#Curated Families:
# AeschNCRlike BowBirk CAMlike CAPCysRich ChiHevWheatin Cystatin Defensin1 Defensin2 EarlyNodulin emp24 HeatShockProt Kunitz LEEDPEED Leginsulin LTP1 LTP2 LTP3 MBP1-1 MBP1-2 NCRs1 NCRs2 NCRs3 NCRs4 NCRs5 NodGRP1 NodGRP2 NodGRP3 NodGRP4 PDP PhaseoleaeNodulin PRP RipeningRelated Thioredoxin


################# SPADA search #################

mkdir -p ~/LSE_pipeline/Stage3/Spada/SPADAprofile/aln/
cd ~/LSE_pipeline/Stage3/Seqs/Curated
#Jalview output alignments must be modified to a correct header format in order to be recognized by SPADA:
for i in AeschNCRlike BowBirk CAMlike CAPCysRich ChiHevWheatin Cystatin Defensin1 Defensin2 EarlyNodulin emp24 HeatShockProt Kunitz LEEDPEED Leginsulin LTP1 LTP2 LTP3 MBP1-1 MBP1-2 NCRs1 NCRs2 NCRs3 NCRs4 NCRs5 NodGRP1 NodGRP2 NodGRP3 NodGRP4 PDP PhaseoleaeNodulin PRP RipeningRelated Thioredoxin; do
  sed -e '1s/CLUSTAL/CLUSTAL 2.1 multiple sequence alignment/' $i.aln > ../../Spada/SPADAprofile/aln/$i.aln
done

#cd ~/software/spada/
perl build_profile.pl --cfg ~/LSE_pipeline/ExampleFiles/Cluster.txt --aln ~/LSE_pipeline/Stage3/Spada/SPADAprofile/aln/ --hmm ~/LSE_pipeline/Stage3/Spada/SPADAprofile/

#Run Spada, final e-value was set to be more stringent for this final search
for i in Tp Mt Pv Gm Ad; do
  perl spada.pl --cfg /home/youngn/truji033/LSE_pipeline/ExampleFiles/Cluster.txt \
    --dir /home/youngn/truji033/LSE_pipeline/Stage3/Spada/$i \
    --hmm /home/youngn/truji033/LSE_pipeline/Stage3/Spada/SPADAprofile/ \
    --fas /home/youngn/truji033/LSE_pipeline/Setup/$i/$i.fa --org $i --evalue 0.1
done

rm -rf ~/LSE_pipeline/Stage3/Spada/*/11_motif_mining ~/LSE_pipeline/Stage3/Spada/*/21_model_prediction

########### Expression of Spada Hits ###########

mkdir -p ~/LSE_pipeline/Stage3/Spada/Expression
cd ~/LSE_pipeline/Stage3/Spada/Expression

#Copy files to run Cufflinks on
#Find tissue expression for the SPADA hits
for j in Tp Mt Pv Gm Ad; do
  mkdir -p $j; cd $j
  cp -s ../../$j/31_model_evaluation/61_final.g* .
  for i in Flower Leaf Root Nodules; do 
    cufflinks -u --min-intron-length 20 --max-intron-length 2000 -G 61_final.gff -p 16 -o $i \
       --max-bundle-frags 50000000 ~/LSE_pipeline/Setup/$j/$i/accepted_hits.sorted.bam 
  done; cd ..
done


####### Tissue expression for Spada hits #######

#Unite Spada IDs to fpkm values for each tissue
for i in Gm Mt Ad Tp Pv; do for j in Flower Leaf Root Nodules; do
awk -F "\t" 'NR==FNR {FPKM[$1]=$10; next} {if (FNR==1) {next}; gsub("*", "", $19); h[$17]+=1 ; \
       print $17"-"h[$17],$19,$2,FPKM[$2]}' $i/$j/genes.fpkm_tracking $i/61_final.gtb > $i/$j-summary.txt
done; done

for i in Gm Mt Ad Tp Pv; do 
  cd $i
  let "num = 1"; let "next = 2"
  for j in Flower Leaf Root Nodules; do
    awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}{print $1}' $j-summary.txt > summary.txt.ALL-1
    awk -v OFS="\t" 'NR==FNR {h[$1]=$4; next} {print $0,h[$1]}' $j-summary.txt summary.txt.ALL-$num > summary.txt.ALL-$next
    let "num = num+1"; let "next = next+1"
  done
  mv summary.txt.ALL-$num ../$i-summary.txt.ALL; rm summary.txt.ALL-*
  cd ..
done


########### Final Stats and Heatmaps ###########

#Create separate fasta and expression files for each gene family
cd ..; mkdir -p Rfigs
for j in AeschNCRlike BowBirk CAMlike CAPCysRich ChiHevWheatin Cystatin Defensin EarlyNodulin emp24 HeatShockProt Kunitz LEEDPEED Leginsulin LTP MBP1 NCRs PDP PhaseoleaeNodulin PRP RipeningRelated Thioredoxin; do
  awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}' > Rfigs/$j.tbl
  for i in Gm Mt Ad Tp Pv; do 
    awk -F "\t" -v Sp=$i -v Grp=$j '{if ($1 ~ Grp) {h[Grp]+=1; print Sp"|"Grp"-"h[Grp]"\t"$2"\t"$3"\t"$4"\t"$5}}' Expression/$i-summary.txt.ALL >> Rfigs/$j.tbl
    awk -F "\t" -v Sp=$i -v Grp=$j '{if (FNR==1) {next}; gsub("*", "", $19);  \
      if ($17 ~ Grp ) {h[Grp]+=1; print ">"Sp"|"Grp"-"h[Grp]"\n" $19}}' Expression/$i/61_final.gtb >> Rfigs/$j.fa
done; done

#Combine GRP groups into two types:
awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}' > Rfigs/GRP1.tbl
awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule"}' > Rfigs/GRP2.tbl
for j in NodGRP1 NodGRP2 NodGRP3; do
  for i in Gm Mt Ad Tp Pv; do 
    awk -F "\t" -v Sp=$i -v Grp=$j '{if ($1 ~ Grp) {h[Grp]+=1; print Sp"|"Grp"-"h[Grp]"\t"$2"\t"$3"\t"$4"\t"$5}}' Expression/$i-summary.txt.ALL >> Rfigs/GRP1.tbl
    awk -F "\t" -v Sp=$i -v Grp=$j '{if (FNR==1) {next}; gsub("*", "", $19); if ($17 ~ Grp ) {h[Grp]+=1; print ">"Sp"|"Grp"-"h[Grp]"\n" $19}}' Expression/$i/61_final.gtb >> Rfigs/GRP1.fa
done; done
for j in NodGRP4; do
  for i in Gm Mt Ad Tp Pv; do 
    awk -F "\t" -v Sp=$i -v Grp=$j '{if ($1 ~ Grp) {h[Grp]+=1; print Sp"|GRP2-"h[Grp]"\t"$2"\t"$3"\t"$4"\t"$5}}' Expression/$i-summary.txt.ALL >> Rfigs/GRP2.tbl
    awk -F "\t" -v Sp=$i -v Grp=$j '{if (FNR==1) {next}; gsub("*", "", $19); if ($17 ~ Grp ) {h[Grp]+=1; print ">"Sp"|GRP2-"h[Grp]"\n" $19}}' Expression/$i/61_final.gtb >> Rfigs/GRP2.fa
done; done

#Align sequences, get dendrogram for figures
for j in AeschNCRlike BowBirk CAMlike CAPCysRich ChiHevWheatin Cystatin Defensin EarlyNodulin emp24 HeatShockProt Kunitz LEEDPEED Leginsulin LTP MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin PRP RipeningRelated Thioredoxin; do
  clustalw2 -infile=Rfigs/$j.fa -type=protein -outfile=Rfigs/$j.aln
done

#Get table with nodule specificity for heatmaps
cd Rfigs
for j in AeschNCRlike BowBirk CAMlike CAPCysRich ChiHevWheatin Cystatin Defensin EarlyNodulin emp24 HeatShockProt Kunitz LEEDPEED Leginsulin LTP MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin PRP RipeningRelated Thioredoxin; do
  awk -v OFS="\t" 'BEGIN{print "ID\tFlower\tLeaf\tRoot\tNodule\tSpecificity"}' > $j.tbl2  
  awk -F "\t" -v Grp=$i '{if ($1!~"ID") {Spec="F"; avgnon=($2+$3+$4)/3; relexp=$5/(avgnon+0.1)
    if ($5>10 && relexp > 4) {Spec="T"}
    print $0"\t"Spec}}' $j.tbl >> $j.tbl2
done

#Combine expression and phylogenetic analyses: Heatmap.R

#Based on the heatmaps-expression analysis, the following showed signs of LSE:
#AeschNCRlike BowBirk CAMlike CAPCysRich Cystatin LEEDPEED Leginsulin MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin

#The following groups did not show signs of LSE:
#ChiHevWheatin Defensin EarlyNodulin emp24 HeatShockProt Kunitz LTP PRP RipeningRelated Thioredoxin


############ Get stats for Table 1 #############  

#Total members per group in each species
for j in AeschNCRlike BowBirk CAMlike CAPCysRich Cystatin LEEDPEED Leginsulin MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin; do
  for i in Gm Pv Mt Tp Ad; do 
  echo $j $i
  grep -c $i $j.tbl
done;done

#Nodule-specific members per group in each species
for j in AeschNCRlike BowBirk CAMlike CAPCysRich Cystatin LEEDPEED Leginsulin MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin; do
  awk '$6=="T"{print $1}' $j.tbl2 >> NodSpecList.txt
done

for j in AeschNCRlike BowBirk CAMlike CAPCysRich Cystatin LEEDPEED Leginsulin MBP1 NCRs GRP1 GRP2 PDP PhaseoleaeNodulin; do
  for i in Gm Pv Mt Tp Ad; do 
  echo $j $i
  grep -c $i"|"$j NodSpecList.txt
done; done




