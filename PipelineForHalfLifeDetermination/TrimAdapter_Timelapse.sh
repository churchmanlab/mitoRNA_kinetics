#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

module load gcc/6.2.0 
module load python/3.7.4 
module load cutadapt/2.5

LibName=$1
fastqName=$2
type=$3 # Truseq Takara

echo $LibName

########## CLEAN READS, TRIM ADAPTERS ###########
if [ "${type}" = "Takara" ]
then
infastq1=${fastqName}_R1.fastq.gz
infastq2=${fastqName}_R2.fastq.gz
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
Adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
fivePtrim="3"
elif [ "${type}" = "Truseq" ]
then # Switch because strandedness is reversed
infastq1=${fastqName}_R2_001.fastq.gz
infastq2=${fastqName}_R1_001.fastq.gz
adapter="CTGTCTCTTATACACATCT"
Adapter="CTGTCTCTTATACACATCT"
fivePtrim="0"
fi

# Cutadapt to trim adapters, allowing 20% error rate since read quality is low at 3' end require minimum length 20 (-m 20), trim 3 nt from 5' end (-u 3), trim Gs that result from sequence runoff [no color, black] (--nextseq-trim=5), max 0 Ns allowed in sequence (max-n 0)
cutadapt -a $adapter -A $Adapter -e .2 -m 20 -u $fivePtrim --nextseq-trim=5 --max-n 0 -j 4 -o ${LibName}_R1_noQtrimmed.fastq -p ${LibName}_R2_noQtrimmed.fastq $infastq1 $infastq2

# Quality at 3' ends is quite low, but adding this quality filter to the above command causes quality trimming before adaptor trimming. In that case about 1,000,000 more reads are discarded, not passing min length filter. Doing it this way, adapter is found in about 3% MORE of the reads, and ~3x fewer bases end up getting quality trimmed
cutadapt -a $adapter -A $Adapter -q 28 -m 20 -j 4 -o ${LibName}_R1_trimmed.fastq -p ${LibName}_R2_trimmed.fastq ${LibName}_R1_noQtrimmed.fastq ${LibName}_R2_noQtrimmed.fastq

# And finally remove 3nt from 3' end of Read2 only to match the 3 trimmed from the 5' end of Read1 (without doing this R2 sometimes extends beyond R1 and causes the pair to not map). Have to take out the min length requirement so no reads get removed and mess up order
cutadapt -j 4 -u -${fivePtrim} -o ${LibName}_R2_trim_trimmed.fastq ${LibName}_R2_trimmed.fastq

if [ "${type}" = "Takara" ]
then
# Also 5 nt from 5' end of R2 to remedy mismatches with random hexamer priming
cutadapt -j 4 -u 5 -o ${LibName}_R2_trim2_trimmed.fastq ${LibName}_R2_trim_trimmed.fastq
# And 5 nt from 3' end of R1 to match
cutadapt -j 4 -u -5 -o ${LibName}_R1_trim_trimmed.fastq ${LibName}_R1_trimmed.fastq

# Remove reads determined to mismap because of concatemerization during lib prep
# RNR2 @1900
../../CytoRP_SMD/Scripts/bbmap/bbduk.sh in=${LibName}_R1_trim_trimmed.fastq in2=${LibName}_R2_trim2_trimmed.fastq hdist=1 k=45 literal=TAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCCCCGAAACC out=${LibName}_R1_trim_trimmed_e.fastq out2=${LibName}_R2_trim2_trimmed_e.fastq
fi



######### Clean up after running on ALL libraries ###########
# rm *_R*_noQtrimmed.fastq
# rm *_R*_trimmed.fastq
 