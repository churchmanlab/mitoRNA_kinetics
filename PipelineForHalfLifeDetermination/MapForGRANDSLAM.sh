#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###


LibName=$1
MapMethod=$2
data=$3 # Hela K562 mouse
seqmeth=$4 # single paired
map_rRNA=$5

echo $LibName

##### To make STAR index: 

# mkdir STAR-index_ensGRCh38_h_MTsnp_fuse_149
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_ensGRCh38_h_MTsnp_fuse_149 --genomeFastaFiles ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta --sjdbGTFfile ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf --sjdbOverhang 149"

# Updated 2/2021 for SNP-corrected genome and refined mito transcript ends
# mkdir STAR-index_Hela_hg38
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_Hela_hg38 --genomeFastaFiles Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta --sjdbGTFfile Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf --sjdbOverhang 149"

# mkdir STAR-index_K562_hg38_dm6
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_K562_hg38_dm6 --genomeFastaFiles K562_ensGRCh38_dm6_ercc_cat.fasta --sjdbGTFfile K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf --sjdbOverhang 149"

# mkdir STAR-index_mouseNIH3T3_mm10_dm6
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_mouseNIH3T3_mm10_dm6 --genomeFastaFiles mouseNIH3T3_mm10_dm6_ercc_cat.fasta --sjdbGTFfile mouseNIH3T3_mm10_MTmod_dm6_ercc_cat.gtf --sjdbOverhang 149"

# Added 7/2022 for LRPPRCKO experiment - SNP-masked from the Hela genome above
# mkdir STAR-index_HEK293T_hg38
# sbatch -p short -t 0-3:00 -c 6 --mem=100G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR-index_HEK293T_hg38 --genomeFastaFiles HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta --sjdbGTFfile HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf --sjdbOverhang 149"


########## define inputs ###########
if [ "${seqmeth}" = "paired" ]
then

infastqR1="${LibName}_R1_trim_trimmed_e.fastq"
infastqR2="${LibName}_R2_trim2_trimmed_e.fastq"

# infastqR1="${LibName}_R1_onlyQtrimmed.fastq"
# infastqR2="${LibName}_R2_onlyQtrimmed.fastq"

# infastqR1="${LibName}_R1_trimmed.fastq"
# infastqR2="${LibName}_R2_trim_trimmed.fastq"

else 
infastqR1="${LibName}_trimmed.fastq"
infastqR2=""
fi

if [ "${data}" = "Hela" ]
then
genome="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/STAR-index_Hela_hg38"
elif [ "${data}" = "K562" ]
then
genome="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/STAR-index_K562_hg38_dm6"
elif [ "${data}" = "mouse" ]
then
genome="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/STAR-index_mouseNIH3T3_mm10_dm6"
elif [ "${data}" = "HEK" ]
then
genome="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/STAR-index_HEK293T_hg38"
fi

########## MAP to rRNA to get estimate of counts (just use one side of pair) ###########
# Optional, uncomment to get ncRNA counts
if [ "${map_rRNA}" = "yes" ]
then
# # Cyto SSU rRNA
bowtie -v 3 -5 2 -3 2 -S /n/groups/churchman/hMitoRP/bowtieindex/rRNA_SSU $infastqR1 ${LibName}_Nuc_SSUrRNA_almts.sam
rm ${LibName}_Nuc_SSUrRNA_almts.sam
# 
# # Cyto LSU rRNA
bowtie -v 3 -5 2 -3 2 -S /n/groups/churchman/hMitoRP/bowtieindex/rRNA_LSU $infastqR1 ${LibName}_Nuc_LSUrRNA_almts.sam
rm ${LibName}_Nuc_LSUrRNA_almts.sam
# 
# # Mito rRNA
bowtie -v 3 -5 2 -3 2 -S /n/groups/churchman/hMitoRP/bowtieindex/hMito_rRNA $infastqR1 ${LibName}_hMito_rRNA_almts.sam
rm ${LibName}_hMito_rRNA_almts.sam
fi

########## MAP ###########
########## note: *trimmed_e.fastq files have known misaligning sequences removed

if [ "${MapMethod}" = "MMinformed3" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_MMinformed3_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 90 --outFilterMismatchNoverLmax 0.3 --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "MMinformed4" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_MMinformed4_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 90 --outFilterMismatchNoverLmax 0.3 --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "t5MMinformed4" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_t5MMinformed4_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 90 --outFilterMismatchNoverLmax 0.3 --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --seedSearchStartLmax 5 --outReadsUnmapped Fastx


elif [ "${MapMethod}" = "stdPlus" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_stdPlus_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 30 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "stdPlus2" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_stdPlus2_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 30 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "MMinformed5" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_MMinformed5_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 15 --outFilterMismatchNoverLmax 0.09 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "t5MMinformed5" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_t5MMinformed5_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 15 --outFilterMismatchNoverLmax 0.09 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "MTMMinformed6" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_MTMMinformed6_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "t5MTMMinformed6" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_t5MTMMinformed6_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

elif [ "${MapMethod}" = "t5MTMMinformed6noTrim" ]
then
	STAR --runThreadN 4 --genomeDir $genome --readFilesIn $infastqR1 $infastqR2 --outFileNamePrefix ${LibName}_t5MTMMinformed6noTrim_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 100 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMultimapScoreRange 0 --seedSearchStartLmax 5 --outReadsUnmapped Fastx

fi


if [ "${data}" = "Hela" ]
then
# Remove ERCC spike ins
samtools view -@ 3 -h ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam | awk '{if ($3 != "ERCC-00048" && $3 != "ERCC-00136") print $0;}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam
# # Index
samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam

elif [ "${data}" = "HEK" ]
then
# Remove ERCC spike ins
samtools view -@ 3 -h ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam | awk '{if ($3 != "ERCC-00048" && $3 != "ERCC-00136") print $0;}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam
# # Index
samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam

elif [ "${data}" = "K562" ]
then
# Remove ERCC spike ins (named differently) and drosophila reads
samtools view -@ 3 -h ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam | awk '{if ($3 != "e_ERCC48" && $3 !~ /^d_/) print $0;}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_Aligned.sortedByCoord.h_noSpike.bam
# # Index
samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam

elif [ "${data}" = "mouse" ]
then
# Remove ERCC spike ins (named differently) and drosophila reads
samtools view -@ 3 -h ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam | awk '{if ($3 != "e_ERCC48" && $3 !~ /^d_/) print $0;}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_Aligned.sortedByCoord.m_noSpike.bam
# # Index
samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam

fi
