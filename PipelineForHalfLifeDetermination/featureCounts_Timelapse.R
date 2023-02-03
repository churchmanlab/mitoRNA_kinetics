#!/n/app/R/4.0.1/bin/Rscript


library(Rsubread)


# Get arguments
args <- commandArgs(trailingOnly = TRUE)

path <- paste0(args[1], '/')
LibName <- args[2]
MapMethod <- args[3]
suffix <- args[4]
ref <- args[5] # Hela K562 mouse
type <- args[6] # multi unique
seqmeth <- args[7]

InputBAMall <- paste0(path, LibName, '_', MapMethod, '_', suffix)
# InputBAMnoDups <- paste0(path, LibName, '_', MapMethod, '_dupsRemoved.bam')

if (ref == 'Hela') {gtffile = '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf'}
if (ref == 'HEK') {gtffile = '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/HEK293T_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf'}
if (ref == 'K562') {gtffile = '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf'}
# if (ref == 'mouse') {gtffile = '
# '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/     '}

if (seqmeth == 'paired') {op = 'TRUE'}
if (seqmeth == 'single') {op = 'FALSE'}

# See https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts


if (type == 'unique') {
DT = featureCounts(InputBAMall, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = FALSE, allowMultiOverlap = TRUE, largestOverlap = TRUE, strandSpecific = 1, isPairedEnd = op, minFragLength = 25, nthreads = 4)
# DTnoDups = featureCounts(InputBAMnoDups, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = FALSE, allowMultiOverlap = TRUE, largestOverlap = TRUE, strandSpecific = 1, isPairedEnd = TRUE, minFragLength = 25, nthreads = 4)
DT_CDS = featureCounts(InputBAMall, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, allowMultiOverlap = TRUE, minOverlap=22, countMultiMappingReads = FALSE, GTF.featureType='CDS', strandSpecific = 1, minFragLength = 25, isPairedEnd = op, nthreads = 4)

write.table(DT$counts, file=paste0(path, LibName,'_', MapMethod, '_featureCounts_unique.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)
# write.table(DTnoDups$counts, file=paste0(path, LibName,'_', MapMethod, '_featureCounts_unique_noDups.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)
write.table(DT_CDS$counts, file=paste0(path, LibName,'_', MapMethod,'_featureCounts_CDSunique.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)

}


if (type == 'multi') {
DT = featureCounts(InputBAMall, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = TRUE, allowMultiOverlap = TRUE, largestOverlap = TRUE, strandSpecific = 1, isPairedEnd = op, minFragLength = 25, nthreads = 4)
# DTnoDups = featureCounts(InputBAMnoDups, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = TRUE, allowMultiOverlap = TRUE, largestOverlap = TRUE, strandSpecific = 1, isPairedEnd = TRUE, minFragLength = 25, nthreads = 4)
DT_CDS = featureCounts(InputBAMall, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, allowMultiOverlap = TRUE, minOverlap=22, countMultiMappingReads = TRUE, GTF.featureType='CDS', strandSpecific = 1, minFragLength = 25, isPairedEnd = op, nthreads = 4)

write.table(DT$counts, file=paste0(path, LibName,'_', MapMethod,'_featureCounts_multi.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)
# write.table(DTnoDups$counts, file=paste0(path, LibName,'_', MapMethod,'_featureCounts_multi_noDups.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)
write.table(DT_CDS$counts, file=paste0(path, LibName,'_', MapMethod,'_featureCounts_CDSmulti.txt'), sep=("\t"), col.names = FALSE, quote=FALSE)

}


# To get feature lengths, for RPK. Just need to do this once
write.table(DT$annotation$Length, file='featureCounts_Length.txt', col.names=FALSE, sep=("\t"), quote=FALSE)
write.table(DT_CDS$annotation$Length, file='featureCounts_CDS_Length.txt', col.names=FALSE, sep=("\t"), quote=FALSE)




