#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com
#SBATCH -e logs/CallSNPs_ModifyGenome_%j.err
#SBATCH -o logs/CallSNPs_ModifyGenome_%j.log


### Written 2020 by M. Couvillion
### This script will call snps with bcftools using one or more bam files, separate the resulting .vcf file into INDELs and non-INDELs, and modify the reference fasta based on the non-INDEL snvs. A custom python script is then used to further split the vcf into loci that have > 75% of coverage as a specific alternate base (these will be changed in the reference to the alternate base), and those that have < 75 % an alternate base (these will be changed to N in the reference). This first reference modification is done with gatk FastaAlternateReferenceMaker. Next a pass filter column is added to the INDEL vcf with another custom python script and another round of reference modification is done from that using rf2m which modifies the matching gtf annotation file in parallel. Finally, matching fasta index, dictionary, and chrom.sizes files are made.


### USE
# Modify this script, then:
# sbatch ../Scripts/CallSNPs_ModifyGenome.sh 



# Load required module
module load picard/2.8.0
module load samtools/1.9
module load bcftools/1.9
module load gatk/4.0.0.0

# helpful tips: https://www.biostars.org/p/327558/

# Set naming paramters
Celline="HEK293T" # Hela K562 mouseNIH3T3 HEK293T
MapMethod="t5MMinformed5" #t5MMinformed5 MMinformed5

# Set path to sequence directory
path="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/"

# Set paths to fasta
infasta='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta' # Use this for new human genomes to match gtf (see below)
# /n/groups/churchman/bms36/genomes/mm10_dm6_ercc_grandslam/mm10_dm6_ercc_cat.fasta
# ensGRCh38_dm6_ercc_cat.fasta 
# ensGRCh38_h_MT_ncRNAs_allERCC_merge.formatted.fasta

out1fasta=${path}ensGRCh38_h_MT_ncRNAs_allERCC_merge.HEKsnpCor.fasta
# mm10_dm6_ercc_cat.snpCor.fasta
# ensGRCh38_dm6_ercc_cat.snpCor.fasta 
# ensGRCh38_h_MT_ncRNAs_allERCC_merge.snpCor.fasta
out2fasta=${path}${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta
# ${Celline}_mm10_dm6_ercc_cat.fasta
# ${Celline}_ensGRCh38_dm6_ercc_cat.fasta
# ${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta

# Set paths to gtf
ingtf=${path}Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf  # Use this for new human genomes since it has anti genes, 7S, etc
# mm10_MTmod_dm6_ercc_cat.gtf 
# ensGRCh38_fusedMTgenes_dm6_ercc_cat.gtf
# ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes2.gtf
outgtf=${path}${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf
# ${Celline}_mm10_MTmod_dm6_ercc_cat.gtf 
# ${Celline}_ensGRCh38_MTmod_dm6_ercc_cat.gtf  
# ${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf

# Need to make sure there is a fasta reference and dictionary for this file in the same directory, e.g.
# faidx reference:
# samtools faidx $infasta
# fasta dictionary
# java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=$infasta O=${infasta/fasta/dict}

# Set path to input files  ** Make sure these have been mapped to infasta above **
inbam1='/n/groups/churchman/mc348/TimelapseSeq/TL10_LRPPRCKO_HEK_2022_07/TL10_0m_WT_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam2="/n/groups/churchman/mc348/TimelapseSeq/TL10_LRPPRCKO_HEK_2022_07/TL10_0m_LRP_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam"


# inbam1='/n/groups/churchman/mc348/TimelapseSeq/Brendan_mouseNIH3T3/G1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam'

# '/n/groups/churchman/mc348/TimelapseSeq/Brendan_mouseNIH3T3/G1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam' 
# "/n/groups/churchman/mc348/TimelapseSeq/Brendan_K562/J1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam" 

# inbam2='/n/groups/churchman/mc348/TimelapseSeq/Brendan_mouseNIH3T3/H1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam'

# '/n/groups/churchman/mc348/TimelapseSeq/Brendan_mouseNIH3T3/H1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam'
# '/n/groups/churchman/mc348/TimelapseSeq/Brendan_mouseNIH3T3/H1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam'
# "/n/groups/churchman/mc348/TimelapseSeq/Brendan_K562/K1_MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam"
# inbam3="/n/groups/churchman/mc348/TimelapseSeq/Hela_2020_09/TL3_0m_t5MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam"
# inbam4="/n/groups/churchman/mc348/TimelapseSeq/Hela_2020_12/TL4_0m_t5MMinformed5_forSNPcalling_Aligned.sortedByCoord.out.bam"

# Call snps 
# This only calls loci with >1 read
# mpileup generates genotype likelihoods at each genomic position with coverage
bcftools mpileup -Ou -f $infasta $inbam1 $inbam2 | bcftools call -mv -Ov -o ${Celline}_${MapMethod}_snpcalls.vcf


# make separate files for INDELs and non-INDELS
grep -v 'INDEL' ${Celline}_${MapMethod}_snpcalls.vcf > ${Celline}_${MapMethod}_snpcalls.noINDEL.vcf
grep -e '^#' -e 'INDEL' ${Celline}_${MapMethod}_snpcalls.vcf > ${Celline}_${MapMethod}_snpcalls.INDELonly.vcf
# make separate files for those that get the alternate allele and those that get N (must have at least 5 reads to get alternate allele)
python /n/groups/churchman/mc348/TimelapseSeq/Scripts/ParseAndFilterVCF.py ${Celline}_${MapMethod}_snpcalls.noINDEL.vcf 0.75
# 
# index noINDEL vcfs for next step
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk IndexFeatureFile -I ${Celline}_${MapMethod}_snpcalls.noINDEL.AFatleast0.75.vcf
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk IndexFeatureFile -I ${Celline}_${MapMethod}_snpcalls.noINDEL.AFbelow0.75.vcf

# First fasta modification without considering INDELs (takes ~30min) 
# This uses gatk to modify fasta with Ns for SNPs with ambiguous alternate allele and the alternate base for SNPs where it's unambiguous
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk FastaAlternateReferenceMaker -R $infasta -O $out1fasta -V ${Celline}_${MapMethod}_snpcalls.noINDEL.AFatleast0.75.vcf --snp-mask ${Celline}_${MapMethod}_snpcalls.noINDEL.AFbelow0.75.vcf --snp-mask-priority true # should this be false?

# Have to delete the number added in front of each fasta header, e.g. 
# >1 1:1-249250621
# >2 10:1-135534747
# >3 h_MT:1-16023
sed -i 's/>.* />/' $out1fasta
sed -i 's/:/ /' $out1fasta

# Remove non-matched index and dictionary files 
rm ${out1fasta/fasta/dict}
rm ${out1fasta/fasta/fasta.fai}

# Next use rf2m to deal with INDELs - this makes a new fasta with matched gtf since coordinates will change when INDELs are allowed

# Need to add PASS to FILTER column for rf2m to read line (number is quality filter cutoff for PASS, output filename ends in '.QUALfilt.vcf')
python /n/groups/churchman/mc348/TimelapseSeq/Scripts/ParseAndPASSlabelVCF.py ${Celline}_${MapMethod}_snpcalls.INDELonly.vcf 10

# Using rf2m to make a new fasta and gtf together, ***including indels***
module load perl/5.30.0

/n/groups/churchman/mc348/programs/rf2m-master/format_vcf.pl --vcf ${Celline}_${MapMethod}_snpcalls.INDELonly.QUALfilt10.vcf

/n/groups/churchman/mc348/programs/rf2m-master/genome_creator.pl --genome $out1fasta --outGenome $out2fasta      ### Had to modify the genome_creator.pl script to read gencode fasta (originally written for ensembl fasta:  $header =~ m/^(.+) dna/; >  $header =~ m/^(.+) 1-/;

/n/groups/churchman/mc348/programs/rf2m-master/gtf_creator.pl --gtf $ingtf --outGTF $outgtf



# Make new fasta index and dictionary
# Make fasta index
samtools faidx $out2fasta
# Make sequence dictionary
java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=$out2fasta O=${out2fasta/fasta/dict}

# Make new chrom.sizes for modified fasta:
cut -f1,2 ${out2fasta}.fai > ${Celline}.chrom.sizes
# To make bed files for separation with samtools, e.g.:
grep '^h_' ${out2fasta}.fai | cut -f1,2 | awk '{print $1 "\t" 0 "\t" $2}' > h_${Celline}_chrNameLength.bed
grep '^d_' ${out2fasta}.fai | cut -f1,2 | awk '{print $1 "\t" 0 "\t" $2}' > d_${Celline}_chrNameLength.bed
grep '^e_' ${out2fasta}.fai | cut -f1,2 | awk '{print $1 "\t" 0 "\t" $2}' > e_${Celline}_chrNameLength.bed