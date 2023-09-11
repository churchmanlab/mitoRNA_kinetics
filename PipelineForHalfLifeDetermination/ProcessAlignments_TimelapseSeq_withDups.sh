#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

###

####SBATCH -c 4
####SBATCH -t 0-24:00
####SBATCH -p medium
####SBATCH --mem=50G
####SBATCH --mail-type=END
####SBATCH --mail-user=mtcouvi@gmail.com

### For small, in silico datasets:
### SBATCH -c 4
### SBATCH -t 0-01:00
### SBATCH -p short
### SBATCH --mem=10G
### SBATCH --mail-type=END
### SBATCH --mail-user=mtcouvi@gmail.com

### Modified from Brendan's processAlignments.sh script

### USE: 			This script will remove secondary alignments and duplicates from the 
###					STAR alignment, remove reads that map to >4 locations and paired
###					reads where the mate did not map properly. The CIGAR string of the
###					remaining reads are then processed so that the mismatches within 
###					each read can be called (soft-clipped bases are removed, and 
###					mismatches are distinguished from matches). Finally, reads are 
###					converted from sam to a pseudo-bed format for quicker downstream
###					analysis.

### REQUIREMENTS:	${libName}_${MapMethod}_${Set}.bam (STAR output, in the same directory)
###  				modifyBed.R (script)
###					path to genome fasta, set below:

libName=$1
MapMethod=$2
Set=$3
data=$4
seqmeth=$5

# Set path to fasta
if [ "${data}" = "Hela" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta'
elif [ "${data}" = "HEK" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta'
elif [ "${data}" = "K562" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta'
elif [ "${data}" = "NIH3T3" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/NIH3T3_mm10_dm6_ercc_cat.fasta'
fi


module load picard/2.8.0




if [ "${seqmeth}" = "paired" ]
then

# sort by read name (required for fixmate)
samtools sort -@ 3 -n -o ${libName}_${MapMethod}_${Set}.nameSort.bam ${libName}_${MapMethod}_${Set}.bam #Aligned.sortedByCoord.subsamp

# Remove reads whose pairs did not make it through the bed filter (this doesn't get rid of the ones that have multiple alignments!)
samtools view -h -@ 3 ${libName}_${MapMethod}_${Set}.nameSort.bam | rev | uniq -D -f16 | rev | samtools view -h -b -S > ${libName}_${MapMethod}_${Set}.nameSort.paired.bam
echo before first unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}.nameSort.bam
echo after first unparied removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}.nameSort.paired.bam

# add fixmate m field (required for markdup)
/n/groups/churchman/bms36/programs/samtools-1.10/bin/samtools fixmate -@ 3 -m ${libName}_${MapMethod}_${Set}.nameSort.paired.bam ${libName}_${MapMethod}_${Set}.nameSort.fm.bam
rm ${libName}_${MapMethod}_${Set}.nameSort.bam
rm ${libName}_${MapMethod}_${Set}.nameSort.paired.bam

# sort by coordinate (required for markdup)
samtools sort -@ 3 -o ${libName}_${MapMethod}_${Set}_2temp.bam ${libName}_${MapMethod}_${Set}.nameSort.fm.bam

else
cp ${libName}_${MapMethod}_${Set}.bam ${libName}_${MapMethod}_${Set}_2temp.bam 

fi

echo Sorted input
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_2temp.bam 



# report read statistics
echo Running samtools flagstat
samtools flagstat -@ 3 ${libName}_${MapMethod}_${Set}_2temp.bam
echo Done

# not primary
samtools view -@ 3 -F 0x100 -o ${libName}_${MapMethod}_${Set}_4temp.bam ${libName}_${MapMethod}_${Set}_2temp.bam
rm ${libName}_${MapMethod}_${Set}_2temp.bam

echo nonprimary removed
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_4temp.bam


# supplementary alignment
samtools view -@ 3 -F 0x800 -o ${libName}_${MapMethod}_${Set}_5temp.bam ${libName}_${MapMethod}_${Set}_4temp.bam
rm ${libName}_${MapMethod}_${Set}_4temp.bam

echo supplementary alignments removed
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_5temp.bam

# reads with unmapped mate
samtools view -@ 3 -F 0x8 -o ${libName}_${MapMethod}_${Set}_6temp.bam ${libName}_${MapMethod}_${Set}_5temp.bam
rm ${libName}_${MapMethod}_${Set}_5temp.bam

echo reads with unmapped mates removed by samtools
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_6temp.bam

# remove reads that map to more than 4 positions, even if it is the primary alignment
samtools view -@ 3 -q 1 -o ${libName}_${MapMethod}_${Set}_filtered_temp.bam ${libName}_${MapMethod}_${Set}_6temp.bam
rm ${libName}_${MapMethod}_${Set}_6temp.bam

echo reads with more than 4 mapping locations removed
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_filtered_temp.bam 

if [ "${seqmeth}" = "paired" ]
then
# # Remove unpaired reads again and should get rid of the rest because they should now be unique (note: out of 10 files with 2-5 mil reads, 2 reads got through. Manually remove them)
# # sort by read name (required for uniq)
samtools sort -@ 3 -n -o ${libName}_${MapMethod}_${Set}_filtered_temp.nameSort.bam ${libName}_${MapMethod}_${Set}_filtered_temp.bam  
# remove unique qnames plus rows that don't have all the attributes (don't have 20 columns)
samtools view -h -@ 3 ${libName}_${MapMethod}_${Set}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | grep '^@\|ms:i:' | samtools view -h -b -S > ${libName}_${MapMethod}_${Set}_filtered_temp.paired.nameSort.bam

echo before second unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_filtered_temp.nameSort.bam
echo after second unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_${Set}_filtered_temp.paired.nameSort.bam

else
cp ${libName}_${MapMethod}_${Set}_filtered_temp.bam ${libName}_${MapMethod}_${Set}_filtered_temp.paired.nameSort.bam

fi

# for removing specific reads too (check logs to make sure no reads escaped, then if so do this)
# samtools view -h -@ 3 ${libName}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | grep -v 'A00794:219:HNH7MDRXX:2:2169:10628:14043\|A00794:219:HNH7MDRXX:2:2169:10628:23281' | samtools view -h -b -S > ${libName}_filtered_temp.paired.nameSort.bam

# sort by coordinate (required for biostar84452.jar)
samtools sort -@ 3 -o ${libName}_${MapMethod}_${Set}_filtered_temp.paired.bam ${libName}_${MapMethod}_${Set}_filtered_temp.paired.nameSort.bam


# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${libName}_${MapMethod}_${Set}_filtered_temp.sam ${libName}_${MapMethod}_${Set}_filtered_temp.paired.bam


# remove soft clipped bases from reads
echo removing soft-clipped bases
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${libName}_${MapMethod}_${Set}_filtered.noSoft_temp.sam ${libName}_${MapMethod}_${Set}_filtered_temp.sam
echo done removing soft-clipped bases

rm ${libName}_${MapMethod}_${Set}_filtered_temp.sam

############ Before running samfixcigar, need to create fasta index and sequence dictionary (just do this once)
############ Load picard: module load picard/2.8.0
############ First, fix fasta so that all lines are the same length (need extra memory for this, I used 50G): java -jar $PICARD/picard-2.8.0.jar NormalizeFasta I=GRCh38_ncRNAs_ERCC_merge.fasta O=GRCh38_ncRNAs_ERCC_merge_linelengthfixed.fasta LINE_LENGTH=60
############ Next, make fasta index: samtools faidx ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta
############ Make picard sequence dictionary: java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta O=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.dict
# converting cigar string to discriminate between mismatches and matches
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/samfixcigar.jar -r $fastapath ${libName}_${MapMethod}_${Set}_filtered.noSoft_temp.sam > ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_${MapMethod}_${Set}_filtered.noSoft_temp.sam

# get rid of sam header
samtools view -@ 3 -O SAM -o ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar.noHead_temp.sam ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.sam

# convert to bam
samtools view -@ 3 -b -o ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.bam ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.sam

# convert to modified bed format
bedtools bamtobed -cigar -i ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.bam > ${libName}_${MapMethod}_${Set}_reads.bed
# rm ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar_temp.bam

echo bed lines
wc -l ${libName}_${MapMethod}_${Set}_reads.bed

# run modifyBed.R script to edit bed file so MM analysis can be run on it *Need more than 25G memory for this step
~/R-3.5.1/library/modifyBed.R ${libName}_${MapMethod}_${Set}_reads.bed ${libName}_${MapMethod}_${Set}_filtered.noSoft.fixCigar.noHead_temp.sam
rm ${libName}_${MapMethod}_${Set}_reads.bed

mkdir ${libName}_${Set}_${MapMethod}_withDups

# split up file to run MM analysis
split --lines=10000 ${libName}_${MapMethod}_${Set}_reads.bedwithC.bed ${libName}_${Set}_${MapMethod}_withDups/tmp.




