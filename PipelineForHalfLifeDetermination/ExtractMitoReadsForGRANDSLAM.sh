#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

module load samtools/1.15.1 # Need to load this for the merge command to work properly

LibName=$1
MapMethod=$2
seqmeth=$3
Set=$4
ref=$5

echo $LibName

if [ "${ref}" = "mouse" ]
then
mito="m_MT"
else
mito="h_MT"
fi


# Index
# samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam

if [ "${Set}" = "MT" ]
then
# Extract mito reads (do once for downstream of rRNA and once for upstream, then combine. Added 9/2/22 to include 7S for mitoriboIP libraries when using downsampled rRNA for GS
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam ${mito}:3305 -o ${LibName}_${MapMethod}_MTdownrRNA.bam
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam ${mito}:1-647 -o ${LibName}_${MapMethod}_MTuprRNA.bam
# merge
samtools merge -f ${LibName}_${MapMethod}_MTnorRNA.bam  ${LibName}_${MapMethod}_MTuprRNA.bam ${LibName}_${MapMethod}_MTdownrRNA.bam


samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam ${mito} -o ${LibName}_${MapMethod}_MT.bam

# # # Index
samtools index -@ 3 ${LibName}_${MapMethod}_MTnorRNA.bam
samtools index -@ 3 ${LibName}_${MapMethod}_MT.bam
elif [ "${Set}" = "Nuc" ]
then
# Remove mito mappers for nuc only
samtools view -@ 3 -h ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam | awk '{if ($3 != "h_MT" && $3 != "m_MT") print $0;}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_Nuc.bam
# # # Index
samtools index -@ 3 ${LibName}_${MapMethod}_Nuc.bam
fi

# Optional for viewing bams:
# Downsample
# samtools view -@ 3 -bs 42.3 ${LibName}_${MapMethod}_MT.bam > ${LibName}_${MapMethod}_MT.ds.bam
# samtools view -@ 3 -bs 42.2 ${LibName}_${MapMethod}_MTnorRNA.bam > ${LibName}_${MapMethod}_MTnorRNA.ds.bam
# 
# # Split into + and - files
# if [ "${seqmeth}" = "paired" ]
# then
# # Downsampled
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "99" || $2 == "355" || $2 == "147" || $2 == "403") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.P.bam
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "83" || $2 == "339" || $2 == "163" || $2 == "419") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.M.bam
# # MTnorRNA
# samtools view -@ 3 -h ${LibName}_${MapMethod}_MTnorRNA.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "99" || $2 == "355" || $2 == "147" || $2 == "403") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MTnorRNA.ds.P.bam
# samtools view -@ 3 -h ${LibName}_${MapMethod}_MTnorRNA.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "83" || $2 == "339" || $2 == "163" || $2 == "419") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MTnorRNA.ds.M.bam
# # else
# # # Downsampled
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "0" || $2 == "256") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.P.bam
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "16" || $2 == "272") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.M.bam
# # # Not downsampled
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MTnorRNA.bam | awk '{if ($1 ~ /^@/ || $2 == "0" || $2 == "256") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MTnorRNA.P.bam
# # samtools view -@ 3 -h ${LibName}_${MapMethod}_MTnorRNA.bam | awk '{if ($1 ~ /^@/ || $2 == "16" || $2 == "272") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MTnorRNA.M.bam
# # fi
# 
# # # Index
# # samtools index ${LibName}_${MapMethod}_MT.ds.P.bam
# # samtools index ${LibName}_${MapMethod}_MT.ds.M.bam
# samtools index ${LibName}_${MapMethod}_MTnorRNA.ds.P.bam
# samtools index ${LibName}_${MapMethod}_MTnorRNA.ds.M.bam
# fi
