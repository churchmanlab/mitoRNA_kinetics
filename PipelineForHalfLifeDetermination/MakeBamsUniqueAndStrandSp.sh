#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

module load samtools/1.15.1 # Need to load this for the merge command to work properly

LibName=$1
MapMethod=$2




# Pull out region of interest
chr="MT"
region=${chr}:600000-660000
if [ "${chr}" != "MT" ] && [ "${chr}" != "MTnorRNA" ]
then
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam $region -o ${LibName}_${MapMethod}_${chr}.bam
fi

# Make file with unique mappers
samtools view -@ 3 -h ${LibName}_${MapMethod}_${chr}.bam | awk '{if ($1 ~ /^@/ || $12 == "NH:i:1") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.uniq.bam

# Separate + and -
samtools view -@ 3 -h ${LibName}_${MapMethod}_${chr}.bam | awk '{if ($1 ~ /^@/ || $2 == "99" || $2 == "355" || $2 == "147" || $2 == "403") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.P.bam
samtools view -@ 3 -h ${LibName}_${MapMethod}_${chr}.bam | awk '{if ($1 ~ /^@/ || $2 == "83" || $2 == "339" || $2 == "163" || $2 == "419") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.M.bam

# Separate + and -
samtools view -@ 3 -h ${LibName}_${MapMethod}_${chr}.uniq.bam | awk '{if ($1 ~ /^@/ || $2 == "99" || $2 == "355" || $2 == "147" || $2 == "403") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.uniq.P.bam
samtools view -@ 3 -h ${LibName}_${MapMethod}_${chr}.uniq.bam | awk '{if ($1 ~ /^@/ || $2 == "83" || $2 == "339" || $2 == "163" || $2 == "419") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.uniq.M.bam

# Index
samtools index ${LibName}_${MapMethod}_${chr}.P.bam
samtools index ${LibName}_${MapMethod}_${chr}.M.bam
samtools index ${LibName}_${MapMethod}_${chr}.uniq.P.bam
samtools index ${LibName}_${MapMethod}_${chr}.uniq.M.bam


# Exp="TL1"
# MapMethod="t5MMinformed4"
# Libs="0m 240m"
# for lib in $Libs
# do
# sbatch ../Scripts/MakeBamsUniqueAndStrandSp.sh ${Exp}_${lib} $MapMethod
# done


