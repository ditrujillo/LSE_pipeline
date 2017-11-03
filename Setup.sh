#!/bin/bash -l

#########################################################################
################################# Setup ################################# 
#########################################################################

#Go to setup folder

cd ~/LSE_pipeline/Setup/
mkdir -p Mt Gm Tp Pv

#Prior to running, the following scripts require manual downloading of genomes from Phytozome
#Connect to JGI Phytozome 
#https://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome#

#Download to ~/LSE_pipeline/Setup/Mt: Mtruncatula_285_Mt4.0.fa.gz 
#Download to ~/LSE_pipeline/Setup/Gm: Gmax_275_v2.0.fa.gz
#Download to ~/LSE_pipeline/Setup/Tp: Tpratense_385_v2.fa.gz
#Download to ~/LSE_pipeline/Setup/Pv: Pvulgaris_442_v2.0.fa.gz
#Then run:
./Mt-Expression.sh
./Gm-Expression.sh
./Tp-Expression.sh
./Pv-Expression.sh

#Arachis duranensis genome Aradu_v1.0_20140908.fa is downloaded from peanutbase.org in the following script:
./Ad-Expression.sh


