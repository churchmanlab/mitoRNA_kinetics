library(data.table)

path <- '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/'
# Input Files
ENStoGeneNameFile <- 'HEK293T_ensGRCh38_ENStoGeneName.txt' # mouseNIH3T3_mm10_ENStoGeneName.txt
exonFile <- 'HEK293T_ensGRCh38_ExonLengths.txt'
# Output Files
maxLengthFile <- 'HEK293T_ensGRCh38_maxLengths.txt'
ENStoGeneNameMaxLength <- 'HEK293T_ensGRCh38_ENStoGeneNameMaxLength.txt'

DTexon <- data.table(read.table(paste0(path, exonFile), sep = '\t', header=TRUE))

# Sum reads over each txpt * keeping GeneName will cause duplications
TxptAggDT <- DTexon[,.(GeneName=GeneName, exonLength.Sum=sum(exonLength)),by=txptID]
# Keep unique txptID lines
TxptAggDT <- TxptAggDT[!duplicated(TxptAggDT$txptID)]
# Now choose max exon length by GeneName
MaxLengthDT <- TxptAggDT[,.(maxLength=max(exonLength.Sum)),by=GeneName]

# Read in ENStoGeneNameFile
DTens <- data.table(read.table(paste0(path, ENStoGeneNameFile), sep = '\t', header=TRUE))

# Merge tables
FinalDT = merge(DTens, MaxLengthDT, by='GeneName')

# Write files
write.table(MaxLengthDT, file=paste0(path, maxLengthFile), row.names=FALSE, sep=("\t"), quote=FALSE)
write.table(FinalDT, file=paste0(path, ENStoGeneNameMaxLength), row.names=FALSE, sep=("\t"), quote=FALSE)

# To run on O2:
# Rscript ../Scripts/SumExonLengthsByTxptChooseMax.R 