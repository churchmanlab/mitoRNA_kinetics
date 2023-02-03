#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###


LibName=$1
MapMethod=$2
include=$3 # both RNR2

echo $LibName


# Index
# samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam


if [ "${include}" = "both" ]
then
# Extract mito rRNA reads 
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam h_MT:648-3305 -o ${LibName}_${MapMethod}_MTrRNA.bam

# Downsample them to 5%
samtools view -@ 3 -bs 42.05 ${LibName}_${MapMethod}_MTrRNA.bam > ${LibName}_${MapMethod}_MTrRNA_ds.bam

# Combine with norRNA bam
samtools merge -@ 3 -f ${LibName}_${MapMethod}_MTdsrRNA.bam ${LibName}_${MapMethod}_MTrRNA_ds.bam ${LibName}_${MapMethod}_MTnorRNA.bam

# Index
samtools index -@ 3 ${LibName}_${MapMethod}_MTdsrRNA.bam
fi



if [ "${include}" = "RNR2" ]
then
# Extract mito rRNA reads 
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam h_MT:1640-3305 -o ${LibName}_${MapMethod}_MTrRNR2.bam

# Downsample them to 10%
samtools view -@ 3 -bs 42.1 ${LibName}_${MapMethod}_MTrRNR2.bam > ${LibName}_${MapMethod}_MTrRNR2_ds.bam

# Combine with norRNA bam
samtools merge -@ 3 -f ${LibName}_${MapMethod}_MTdsRNR2.bam ${LibName}_${MapMethod}_MTrRNR2_ds.bam ${LibName}_${MapMethod}_MTnorRNA.bam

# Index
samtools index -@ 3 ${LibName}_${MapMethod}_MTdsRNR2.bam
fi

