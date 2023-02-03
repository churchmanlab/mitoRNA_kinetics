#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

# First run MakeBamsUniqueAndStrandSp.sh
# And after this move to PC and view on IGV or run SplitReadsPercent.r

module load samtools/1.15.1 # Need to load this for the merge command to work properly

LibName=$1
MapMethod=$2




# Pull out region of interest
chr='MT'

infileplus=${LibName}_${MapMethod}_${chr}.P.bam

infileminus=${LibName}_${MapMethod}_${chr}.M.bam


# Header is modified here to only list h_MT

# Separate "blue" and "orange"
samtools view -@ 3 -h $infileplus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/  || ($17 ~ /T/ && $17 !~ /A/) ) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.P.blue.bam
samtools view -@ 3 -h $infileplus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/  || ($17 ~ /A/ && $17 !~ /T/)) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.P.orange.bam
samtools view -@ 3 -h $infileplus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/  || ($17 ~ /A/ && $17 ~ /T/) || ($17 !~ /A/ && $17 !~ /T/ && $1 !~ /^@/) ) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.P.neither.bam
# 
samtools view -@ 3 -h $infileminus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/ || ($17 ~ /T/ && $17 !~ /A/)) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.M.blue.bam
samtools view -@ 3 -h $infileminus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/  || ($17 ~ /A/ && $17 !~ /T/) )print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.M.orange.bam
samtools view -@ 3 -h $infileminus | awk '{if ($1 ~ /^@HD/ || $2 ~ /SN:h_MT/ || $1 ~ /^@PG/ || $1 ~ /^@CO/  || ($17 ~ /A/ && $17 ~ /T/) || ($17 !~ /A/ && $17 !~ /T/ && $1 !~ /^@/) ) print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_${chr}.M.neither.bam
# 
# # Index
samtools index ${LibName}_${MapMethod}_${chr}.P.blue.bam
samtools index ${LibName}_${MapMethod}_${chr}.P.orange.bam
samtools index ${LibName}_${MapMethod}_${chr}.P.neither.bam

samtools index ${LibName}_${MapMethod}_${chr}.M.blue.bam
samtools index ${LibName}_${MapMethod}_${chr}.M.orange.bam
samtools index ${LibName}_${MapMethod}_${chr}.M.neither.bam
# 
# # Make pileups
genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.P.blue.bam > ${LibName}_${MapMethod}_${chr}.P.blue.bedGraph
genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.P.orange.bam > ${LibName}_${MapMethod}_${chr}.P.orange.bedGraph
genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.P.neither.bam > ${LibName}_${MapMethod}_${chr}.P.neither.bedGraph

genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.M.blue.bam > ${LibName}_${MapMethod}_${chr}.M.blue.bedGraph
genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.M.orange.bam > ${LibName}_${MapMethod}_${chr}.M.orange.bedGraph
genomeCoverageBed -bga -trackline -split -ibam ${LibName}_${MapMethod}_${chr}.M.neither.bam > ${LibName}_${MapMethod}_${chr}.M.neither.bedGraph

# Exp="TL3" # TL1
# MapMethod="t5MTMMinformed6" # t5MMinformed4
# Libs="0m 240m"
# for lib in $Libs
# do
# sbatch ../Scripts/BinByConversion.sh ${Exp}_${lib} $MapMethod
# done


