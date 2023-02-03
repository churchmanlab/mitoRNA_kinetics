#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

Exp=$1
Lib=$2
MapMethod=$3
MTstatus=$4
region=$5
count=$6
filt=$7

echo $Lib

module load picard/2.8.0


inputbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_filtered_temp.paired.nameSort.bam
namesortbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_filtered_picardnamesort.bam

if [ "${filt}" = "dist" ] 
then
filterbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_maxTcount${count}.bam
coordsortbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_maxTcount${count}_sort.bam
readlist=MMfrequency/${Exp}_${MTstatus}_${MapMethod}_withDups_frag_${region}_readNames_${Exp}_${Lib}maxTcount${count}.txt

elif [ "${filt}" = "numTC" ]
then
filterbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_GrThan${count}TCperRead.bam
coordsortbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_GrThan${count}TCperRead_sort.bam
readlist=MMfrequency/${Exp}_${MTstatus}_${MapMethod}_withDups_frag_${region}_readNames_${Exp}_${Lib}GrTh${count}TCperRead.txt

elif [ "${filt}" = "numOther" ]
then
filterbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_${count}orLessOtherMM.bam
coordsortbam=${Exp}_${Lib}_${MapMethod}_${MTstatus}_${count}orLessOtherMM_sort.bam
readlist=MMfrequency/${Exp}_${MTstatus}_${MapMethod}_withDups_frag_${region}_readNames_${Exp}_${Lib}_${count}orLessOtherMM.txt

fi


# Sort by readname
java -jar $PICARD/picard-2.8.0.jar SortSam I=$inputbam O=$namesortbam SORT_ORDER=queryname

# Filter by readname
java -jar $PICARD/picard-2.8.0.jar FilterSamReads I=$namesortbam O=$filterbam READ_LIST_FILE=$readlist SORT_ORDER=queryname FILTER=includeReadList

# Sort by coordinate
samtools sort -@ 3 -o $coordsortbam $filterbam

# Index
samtools index -@ 3 $coordsortbam


########## Make bamlist AFTER running above on all samples ##########
# Exp='TL9' 
# MapMethod='t5MTMMinformed6' # MTMMinformed6 t5MTMMinformed6
# MTstatus='MTnorRNA'
# maxcount='29'
# # 
# bglib="0m_tot"
# no4sUlib="no4sU_tot"
# # # # # Copy 0m sample to make 4sU sample for comparison in GS
# cp ${Exp}_${bglib}_${MapMethod}_${MTstatus}_maxcount${maxcount}_sort.bam ${Exp}_${no4sUlib}_${MapMethod}_${MTstatus}_maxcount${maxcount}_sort.bam
# cp ${Exp}_${bglib}_${MapMethod}_${MTstatus}_maxcount${maxcount}_sort.bam.bai ${Exp}_${no4sUlib}_${MapMethod}_${MTstatus}_maxcount${maxcount}_sort.bam.bai
# 
# # Make bamlist
# Experiment="TL9_tot" # TL5_poly
# touch ${Experiment}_${MTstatus}_${MapMethod}_maxcount${maxcount}.bamlist
# # Libs="no4sU_IP 0m_IP 15m_IP 30m_IP 60m_IP"
# Libs="no4sU_tot 0m_tot 15m_tot 30m_tot 60m_tot"
# 
# for lib in $Libs
# do
# echo "${Exp}_${lib}_${MapMethod}_${MTstatus}_maxcount${maxcount}_sort.bam" >> ${Experiment}_${MTstatus}_${MapMethod}_maxcount${maxcount}.bamlist
# done 
# 

# # ########## For quick batch submission ##########
# Exp='TL11' 
# Libs="0m_tot 15m_tot 30m_tot 60m_tot 0m_IP 15m_IP 30m_IP 60m_IP"
# MapMethod='t5MMinformed4' # MTMMinformed6 t5MTMMinformed6
# MTstatus='MTnorRNA'
# region='MTnorRNA'
# count='1'
# for lib in $Libs
# do
# sbatch -e logs/ExtractSelect_${lib}_${count}.err -o logs/ExtractSelect_${lib}_${count}.log ../Scripts/ExtractSelectReadsByNameForGRANDSLAM.sh $Exp $lib $MapMethod $MTstatus $region $count
# done

