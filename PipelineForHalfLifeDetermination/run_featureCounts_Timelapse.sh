#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com
###


LibName=$1
MapMethod=$2
suffix=$3
ref=$4 # Hela K562 mouse HEK
type=$5 # unique multi
seqmeth=$6
path=`pwd`


module load R/4.0.1


echo $LibName


if [ "$ref" = "Hela" ]
then
# Run featureCounts
/n/groups/churchman/mc348/TimelapseSeq/Scripts/featureCounts_Timelapse.R $path $LibName $MapMethod ${suffix} $ref $type $seqmeth

elif [ "$ref" = "K562" ]
then
# Sort and index
samtools sort -@ 3 -o ${LibName}_${MapMethod}_sort.${suffix} ${LibName}_${MapMethod}_${suffix}
samtools index -@ 3 ${LibName}_${MapMethod}_sort.${suffix}
# Run featureCounts
/n/groups/churchman/mc348/TimelapseSeq/Scripts/featureCounts_Timelapse.R $path $LibName $MapMethod sort.${suffix} $ref $type $seqmeth

elif [ "$ref" = "HEK" ]
then
# Run featureCounts
/n/groups/churchman/mc348/TimelapseSeq/Scripts/featureCounts_Timelapse.R $path $LibName $MapMethod ${suffix} $ref $type $seqmeth

fi

