#!/bin/bash

#SBATCH -c 1
#SBATCH -t 360
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH --mail-type=END
#SBATCH --mail-user=email@gmail.com

module load gcc

prefix=$1

echo "Mismatches (Genome nt, Read nt):"
echo "genome_read="$prefix
/n/groups/churchman/mc348/TimelapseSeq/Scripts/countMM_bed.R reads.MM.bed_sort.bed $prefix # reads_MM.fragments.sort.bed



############## combine files ############
# Libs="DMSO 1_NHC 10_NHC 100_NHC 100_Molnu"
# for lib in $Libs
# do
# sed -i '1d' ${Exp}_${lib}_${dirName}/${lib}_mismatches.count
# done
# 
# 
# Exp="SStot1"
# lib=("DMSO" "1_NHC" "10_NHC" "100_NHC" "100_Molnu")
# dirName='Nuc_MMinformed4_withDups'
# 
# paste <(awk 'BEGIN { FS = "=" } {print $1 "\t"$2}' ${Exp}_${lib[0]}_${dirName}/${lib[0]}_mismatches.count) <(awk 'BEGIN { FS = "=" } {print $2}' ${Exp}_${lib[1]}_${dirName}/${lib[1]}_mismatches.count) <(awk 'BEGIN { FS = "=" } {print $2}' ${Exp}_${lib[2]}_${dirName}/${lib[2]}_mismatches.count) <(awk 'BEGIN { FS = "=" } {print $2}' ${Exp}_${lib[3]}_${dirName}/${lib[3]}_mismatches.count) <(awk 'BEGIN { FS = "=" } {print $2}' ${Exp}_${lib[4]}_${dirName}/${lib[4]}_mismatches.count) > ${Exp}_${dirName}_mismatches.count 




