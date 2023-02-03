#!/bin/bash

#SBATCH -c 1
#SBATCH -t 70
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH -o logs/run_AWKforTCperT_%j.out
#SBATCH -e logs/run_AWKforTCperT_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com


### USE: 		Run from one directory up from directories with reads_MM.fragments.sort.bed file 
###				and make MismatchFrequency directory	

### REQUIREMENTS:	MismatchFrequencyPerTplot.R

# Paired-end SAM strand specification: https://biobeat.wordpress.com/2013/04/29/directional-rna-seq-part-1-extract-strand-information-from-sam-file/
# Top strand: 99 = R1+,  147 = R2- (primary), 355 = R1+, 403 = R2- (not primary)
# Bottom strand: 83 = R1-, 163 = R2+, (primary), 339 = R1-, 419 = R2+ (not primary)



LibName=$1
MapMethod=$2
seqmeth=$3
second=$4
Set=$5 # MT MTnorRNA NucTop2000

if [ "${seqmeth}" = "paired" ]
then
infile='reads_MM.fragments.sort.bed'
outname='frag_MMfrequency'
else
infile='reads.MM.bed_sort.bed'
outname='read_MMfrequency'
fi


# for plus strand only: if (($7 == "99" || $7 == "355" || $7 == "0") && $1 != "h_MT"
if [ "${second}" = "CT" ]
then

if [ "${Set}" == "MT" ]
then
# For all mito reads
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tread_length\ttot_mismatches\tTC_mismatches\tCT_mismatches\treadName";} {if ($1 == "h_MT" || $1 == "m_MT") print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$23"\t"$4}' ${LibName}_${MapMethod}/$infile > ${LibName}_${MapMethod}/${outname}_${Set}.txt
fi


if [[ "${Set}" == *"Nuc"* ]]
then
# For all nuc reads
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tread_length\ttot_mismatches\tTC_mismatches\tCT_mismatches\treadName";} {if ($1 != "h_MT" && $1 != "m_MT") print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$23"\t"$4}' ${LibName}_${MapMethod}/$infile > ${LibName}_${MapMethod}/${outname}_${Set}.txt
fi

if [ "${Set}" == "MTnorRNA" ]
then
# # # mito reads, skip rRNAs
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches\tCT_mismatches\treadName";} {if (($1 == "h_MT" || $1 == "m_MT") && ($3 < 648 || $3 > 3305)) print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$23"\t"$4}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_${Set}.txt
fi

fi




if [ "${second}" = "GA" ]
then

if [ "${Set}" == "MT" ]
then
# For all mito reads
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tread_length\ttot_mismatches\tTC_mismatches\tGA_mismatches\treadName";} {if ($1 == "h_MT" || $1 == "m_MT") print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$20"\t"$4}' ${LibName}_${MapMethod}/$infile > ${LibName}_${MapMethod}/${outname}_${Set}.txt
fi



if [[ "${Set}" == *"Nuc"* ]]
then
# For all nuc reads
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tread_length\ttot_mismatches\tTC_mismatches\tGA_mismatches\treadName";} {if ($1 != "h_MT" && $1 != "m_MT") print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$20"\t"$4}' ${LibName}_${MapMethod}/$infile > ${LibName}_${MapMethod}/${outname}_${Set}.txt
fi

if [ "${Set}" == "MTnorRNA" ]
then
# # # mito reads, skip rRNAs
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches\tGA_mismatches\treadName";} {if (($1 == "h_MT" || $1 == "m_MT") && ($3 < 648 || $3 > 3305)) print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$(NF)"\t"$18"\t"$19"\t"$20"\t"$4}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_${Set}.txt
fi

fi

