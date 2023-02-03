#!/bin/bash

#SBATCH -c 4
#SBATCH -e logs/GS_%j.err
#SBATCH -o logs/GS_%j.log
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

######
####SBATCH --x11=batch


ExpName=$1 # TL4_MT TL3_MT TL3_SpikeOnly MUT-TL3_0m-0.1 K562_J_MT TL5_tot_MT TL7_MT TL8_tot_MT TL8_IP_MT TL3_noSpikeNoMito TL8_tot_MTdsrRNA SStot1_Nuc SStot1_MTnorRNA_MDmodC TL9_tot_MT
MapMethod=$2 # t5MTMMinformed6_MT MMinformed4 t5MTMMinformed6 t5MMinformed5
data=$3 # Hela K562 HelaOld HelaCtoT
Set=$4
makeCIT=$5
dirName=$6 # PcModInput PcModErik1 ErikPcTot1

module load gcc/6.2.0
module load java/jdk-1.8u112
module load R/4.0.1
export PATH=$PATH:/n/groups/churchman/mc348/TimelapseSeq/Scripts/GRAND-SLAM_2.0.5f
# export PATH=$PATH:/n/groups/churchman/mc348/TimelapseSeq/Scripts/GRAND-SLAM_2.0.5g-preview


############## First build genome reference #################
########### Notes. This doesn't need to be redone ###########

# Try running everything with ensGRCh38 again, but this time output sorted bam from STAR instead of sorting after

# 1. change name MT to h_MT in fasta and gtf
# sed 's/MT/h_MT/' ensGRCh38_ncRNAs_ERCC_merge.formatted.fasta > ensGRCh38_h_MT_ncRNAs_ERCC_merge.formatted.fasta
# sed 's/^MT/h_MT/' ensGRCh38_ncRNAs_ERCC_merge.gtf > ensGRCh38_h_MT_ncRNAs_ERCC_merge.gtf

# 2. find snps in mito with bcftools 2/18/21: redid this for whole genome and modified gtf to match
# See CallSNPs.sh (use bam file previously mapped to ensGRCh38 with MMinformed3, sorted)

# 3. replace spike-in ERCCs with all ERCCs in fasta and gtf
# Replace the two ERCC transcripts in fasta and gtf with the full set:
# Remove last 37 lines from fasta (these are ERCC) 51663764-37= 51663727
# head -n 51663727 ensGRCh38_h_MT_ncRNAs_ERCC_merge.formatted.fasta > ensGRCh38_h_MT_ncRNAs_ERCCtemp_merge.formatted.fasta
# cat ensGRCh38_h_MT_ncRNAs_ERCCtemp_merge.formatted.fasta ERCC92.fa > ensGRCh38_h_MT_ncRNAs_allERCC_merge.formatted.fasta
# Remove last 2 lines from gtf (ERCC) 2905064-2=2905062
# head -n 2905062 ensGRCh38_h_MT_ncRNAs_ERCC_merge.gtf > ensGRCh38_h_MT_ncRNAs_ERCCtemp_merge.gtf
# cat ensGRCh38_h_MT_ncRNAs_ERCCtemp_merge.gtf ERCC92.gtf > ensGRCh38_h_MT_ncRNAs_allERCC_merge.gtf

# 4. fuse ATP8/6 and ND4L/4 in gtf
# Using vim, see ensGRCh38_h_MTsnp_ncRNAs_ERCC_merge_fusedMTgenes.README.txt

# 5. update mito snps in fasta 2/18/21 Now do this automatically with gatk and rf2m

# 6. make new STAR index
# mkdir STAR-index_ensGRCh38_h_MTsnp_fuse_149
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_ensGRCh38_h_MTsnp_fuse_149 --genomeFastaFiles ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta --sjdbGTFfile ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf --sjdbOverhang 149"
# 2/18/21
# mkdir STAR-index_Hela_hg38
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_Hela_hg38 --genomeFastaFiles Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta --sjdbGTFfile Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf --sjdbOverhang 149"

# 7. make new GS index
# mkdir GRAND_SLAM_ensGRCh38_h_MTsnp_fuse_allERCC_aMito
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf -n ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge_fuse_aMito -D -f GRAND_SLAM_ensGRCh38_h_MTsnp_fuse_allERCC_aMito -nobowtie -nokallisto -nostar"

# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_ncRNAs_allERCC_mitoCirc_merge_fusedMTgenes.gtf -n ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge_fuse_aMito -D -f GRAND_SLAM_ensGRCh38_h_MTsnp_fuse_allERCC_aMito -nobowtie -nokallisto -nostar"

# 10/2020
# Found a deletion in the mito genome the wasn't fixed with SNPs already so add that in and remake indexes
# h_MT 3107	TACNTTC (delete N)
# 1. Update mito snp in fasta 
# 2. Fix line length in fasta
# module load picard/2.8.0
# sbatch -p short -t 0-01:00 --mem=5G --wrap="java -jar $PICARD/picard-2.8.0.jar NormalizeFasta I=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.fasta O=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta LINE_LENGTH=60"
# 3. make new STAR index
# mkdir STAR-index_ensGRCh38_h_MTsnp_fuse_149
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_ensGRCh38_h_MTsnp_fuse_149 --genomeFastaFiles ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta --sjdbGTFfile ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf --sjdbOverhang 149"
# 4. make new GS index
# mkdir GRAND_SLAM_ensGRCh38_h_MTsnp_fuse_allERCC_aMito
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf -n ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge_fuse_aMito -D -f GRAND_SLAM_ensGRCh38_h_MTsnp_fuse_allERCC_aMito -nobowtie -nokallisto -nostar"

# 2/18/21
# mkdir GRAND_SLAM_Hela_ensGRCh38_h_MT
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf -n Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod -D -f GRAND_SLAM_Hela_ensGRCh38_h_MT -nobowtie -nokallisto -nostar"

# 4/10/22
# Make genome where all Cs are changed to Ts
# sed -e '/^>/!s/T/C/g' Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta > Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_CtoT.fasta # in interactive session
# mkdir GRAND_SLAM_Hela_ensGRCh38_h_MT_CtoT
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_CtoT.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf -n Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod_CtoT -D -f GRAND_SLAM_Hela_ensGRCh38_h_MT_CtoT -nobowtie -nokallisto -nostar"


# remade 9/15/22 after adding 7S RNA to gtf
# mkdir GRAND_SLAM_K562_ensGRCh38_dm6 
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf -n K562_ensGRCh38_MTmod_dm6_ercc_cat -D -f GRAND_SLAM_K562_ensGRCh38_dm6 -nobowtie -nokallisto -nostar"

# 7/19/22
# mkdir GRAND_SLAM_HEK_ensGRCh38_h_MT
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf -n HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod -D -f GRAND_SLAM_HEK293T_ensGRCh38_h_MT -nobowtie -nokallisto -nostar"

# 9/15/22 mouse
# mkdir GRAND_SLAM_mouse3T3_mm10_dm6 
# sbatch -p short -t 0-06:00 --mem=80G --wrap="gedi -e IndexGenome -s /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_dm6_ercc_cat.fasta -a /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_MTmod_dm6_ercc_cat.gtf -n mouseNIH3T3_mm10_MTmod_dm6_ercc_cat -D -f GRAND_SLAM_mouse3T3_mm10_dm6 -nobowtie -nokallisto -nostar"

#############################################################
#############################################################
if [ "${data}" = "Hela" ]
then
genomeIndex="Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod"
elif [ "${data}" = "HelaCtoT" ]
then
genomeIndex="Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod_CtoT" # This was removed
elif [ "${data}" = "HelaOld" ]
then
genomeIndex="ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge_fuse_aMito"
elif [ "${data}" = "K562" ]
then
genomeIndex="K562_ensGRCh38_MTmod_dm6_ercc_cat"
elif [ "${data}" = "HEK" ]
then
genomeIndex="HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod"
elif [ "${data}" = "mouse" ]
then
genomeIndex="mouseNIH3T3_mm10_MTmod_dm6_ercc_cat"
fi

# Make CIT files 
if [ "${makeCIT}" = "yes" ] 
then
bamlist2cit -p ${ExpName}_${MapMethod}.bamlist 
fi

if [ "${Set}" = "Nuc" ] 
then
# Run GS ***** Nuc
mkdir GS_${ExpName}_${MapMethod}_lenient_modeAll
gedi -e Slam -allGenes -full -lenient -nthreads 4 -mode All -genomic $genomeIndex -prefix GS_${ExpName}_${MapMethod}_lenient_modeAll/${ExpName}_${MapMethod}_lenient_modeAll -progress -strandness Sense -D -snpConv 0.4 -trim5p 10 -trim3p 5 -reads ${ExpName}_${MapMethod}.bamlist.cit

elif [ "${Set}" = "NucRates" ]
then

# Run GS on modeled rates ***** Nuc
gedi -e Slam -allGenes -full -lenient -nthreads 4 -mode All -genomic $genomeIndex -prefix GS_${ExpName}_${MapMethod}_lenient_modeAll_${dirName}/${ExpName}_${MapMethod}_lenient_modeAll -progress -strandness Sense -D -snpConv 0.4 -trim5p 10 -trim3p 5 -reads ${ExpName}_${MapMethod}.bamlist.cit

elif [ "${Set}" = "MT" ]
then

# Run GS (without -lenient) ***** MT
mkdir GS_${ExpName}_${MapMethod}_modeAll
gedi -e Slam -allGenes -full -nthreads 4 -mode All -genomic $genomeIndex -prefix GS_${ExpName}_${MapMethod}_modeAll/${ExpName}_${MapMethod}_modeAll -progress -strandness Sense -D -snpConv 0.4 -trim5p 10 -trim3p 5 -reads ${ExpName}_${MapMethod}.bamlist.cit

elif [ "${Set}" = "MTrates" ]
then

# Run GS modified rates # delete .tsv and ext.tsv and modify .rates.tsv ****** MT
gedi -e Slam -allGenes -full -nthreads 4 -mode All -genomic $genomeIndex -prefix GS_${ExpName}_${MapMethod}_modeAll_${dirName}/${ExpName}_${MapMethod}_modeAll -progress -strandness Sense -D -snpConv 0.4 -trim5p 10 -trim3p 5 -reads ${ExpName}_${MapMethod}.bamlist.cit

fi

# Usage:
# gedi -e Slam <Options>
# 
# General:
#  -prefix <prefix>                        The prefix used for all output files
#  -nthreads <nthreads>                    The number of threads to use for computations (default: 8)
#  -genomic <genomic>                      The indexed GEDI genome.
#  -reads <reads>                          The mapped reads from the SLAM-seq experiment.
#  -trim5p <trim5p>                        The number bases to trim at the 5' ends of reads (default: 0)
#  -trim3p <trim3p>                        The number bases to trim at the 3' ends of reads (default: 0)
#  -no4sUpattern <no4sUpattern>            The pattern to recognize no4sU conditions (default: no4sU|nos4U)
# 
# SlamDetectSnps:
#  -snpConv <snpConv>                      Conversion rate for SNP calling (default: 0.3)
#  -snppval <snppval>                      Minimal posterior probability of being a conversion (default: 0.001)
# 
# SlamCollectStatistics:
#  -allGenes                               Keep all genes
#  -nUMI                                   Compute nUMIs
#  -mode <mode>                            Read count mode what to do with multi-mapping reads (default: Weight) Choices are : All,Collapse,Divide,Unique,CollapseUnique,Weight
#  -overlap <overlap>                      Overlapping gene mode: What to do for locations that are compatible with more than one gene (default: All) Choices are : All,Collapse,Divide,Unique,CollapseUnique,Weight
#  -introns                                Compute introns
#  -lenient                                Apply lenient transcript compatibility (the two leading and trailing nt may be intronic, reads may extend beyond transcript boundaries); use this for 10x CellRanger mapped reads!
# 
# SlamRates:
#  -plot                                   Use R to produce various plots
#  -minEstimateReads <minEstimateReads>    Specify the minimal number of reads to be used for parameter estimation! (default: 10000)
#  -errlm <errlm>                          Specify the error rate correction linear model (which is applied to the other 11 error rates; expression is evaluated by JS, variables are AC, AG, AT, ... and median) (default: TA*1.434)
#  -errlm2 <errlm2>                        Specify the error rate correction linear model (which is applied to the other 11 error rates; expression is evaluated by JS, variables are AC, AG, AT, ... and median) (default: TA*1.434)
# 
# SlamInfer:
#  -lower <lower>                          The lower credible interval bound (default: 0.05)
#  -upper <upper>                          The upper credible interval bound (default: 0.95)
#  -full                                   Full output (including credible interval, coverages etc.
#  -sparse                                 write sparse matrix format
#  -double                                 Use only doubly sequenced parts for estimation
# 
# SlamTest:
#  -test                                   A test
# 
# Commandline:
#  -progress                               Show progress
#  -D                                      Verbose output of errors
#  -h                                      Show usage
#  -hh                                     Show verbose usage
#  -hhh                                    Show extra verbose usage
#  -dry                                    Dry run
#  -keep                                   Do not remove temp files
