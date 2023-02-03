library(data.table)
library(rlist)

method = 'featureCounts' # star featureCounts
type = 'CDSmulti' # multi CDSmultiMT_multi_noDups

if (method == 'star') {type = 'unique'}

pseudo = 'exclude' # include exclude
ncRNA = 'exclude'
Folder = 'SteadyState2019' #'Hela_2020_09' Brendan_K562_InSilico Brendan_K562 Brendan_mouseNIH3T3 Hela_2021_03_poly Hela_TL4_combined SteadyState2019 Hela_2021_TL6_doubleLabeling Hela_TL3_2020_09 SStot1_COVIDdrugTest_2022_02 SStot2_COVIDdrugTest_2022_04 Brendan_K562_50 Hela_TL2_4sUtest_2020_07 NHC2_HelaTimecourse_2022_06 HEK_TL10_LRPPRC_2022_07 Hela_TL11_500IP_2022_07 Hela_NHC3_500IP_2022_07
Exp = 'SS1' # K562_K TL5_poly SStot2
MapMethod = 'MTMMinformed6' # MMinformed4 t5MTMMinformed6 MTMMinformed6 MMinformed5
genome = 'Hela'



if (method == 'star') {
fileName = paste0(Exp, '_', MapMethod, '_ReadsPerGene')} #'TL3_t5MTMMinformed6_ReadsPerGene' 
if (method == 'featureCounts') {
fileName = paste0(Exp, '_', MapMethod, '_featureCounts_', type)} #'TL3_t5MTMMinformed6_featureCounts_' 
path = paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',Folder,'/RPK/')

DT <- data.table(read.table(paste0(path, fileName, '.txt'), sep = '\t', quote = '', header=TRUE, stringsAsFactors = FALSE))

if (genome == 'Hela') {
GeneConvTable <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Hela_ensGRCh38_ENStoGeneNameMaxLength.txt', header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE))}

if (genome == 'HEK') {
GeneConvTable <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/HEK293T_ensGRCh38_ENStoGeneNameMaxLength.txt', header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE))}

if (genome == 'K562') {
GeneConvTable <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/K562_ensGRCh38_ENStoGeneNameMaxLength.txt', header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE))}

if (genome == 'mouse') {
GeneConvTable <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_ENStoGeneNameMaxLength.txt', header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE))}

GeneConvTableNoPseudo <- GeneConvTable[!grep('pseudogene', gene_biotype)] 
GeneConvTableNoPseudo <- GeneConvTableNoPseudo[!grep('MTRNR2L', GeneName)] 
GeneConvTableNoRNA <- GeneConvTable[!grep('RNA', gene_biotype)] 
GeneConvTableNoPseudoNoRNA <- GeneConvTableNoPseudo[!grep('RNA', gene_biotype)]

# Merge by GeneID(stable), only adding the GeneName and Length columns from the conversion table
if (pseudo == 'include' & ncRNA == 'include'){
mergeDT <- merge(DT, GeneConvTable[, c('GeneID', 'GeneName', 'maxLength')], by='GeneID', all = FALSE)
}
if (pseudo == 'exclude' & ncRNA == 'include'){
mergeDT <- merge(DT, GeneConvTableNoPseudo[, c('GeneID', 'GeneName', 'maxLength')], by='GeneID', all = FALSE)
}
if (pseudo == 'exclude' & ncRNA == 'exclude'){
mergeDT <- merge(DT, GeneConvTableNoPseudoNoRNA[, c('GeneID', 'GeneName', 'maxLength')], by='GeneID', all = FALSE)
}

if (method == 'star') {
SampNum=ncol(mergeDT)-3
samplenames = colnames(mergeDT)[2:(SampNum+1)]
setcolorder(mergeDT, c('GeneName',samplenames,'maxLength'))
setnames(mergeDT, 'maxLength', 'Length')}

if (method == 'featureCounts') {
SampNum=ncol(mergeDT)-4
samplenames = colnames(mergeDT)[2:(SampNum+1)]
setcolorder(mergeDT, c('GeneName',samplenames,'Length'))
mergeDT[, maxLength := NULL]}

# Keep only mito genes, including ncRNAs
# DT <- DT[grep('^MT-', GeneName)]

###### Following are examples of calculating normalization numbers directly from this file
# x=samplecolumns
# for (i in x){
# si_nuc = (sum(mDT[[i]])-sum(mDT[grep('^mt-', GeneName)][[i]]))/(sum(mDT[[i]]) + sum(hDT[[i]])) *100
# si_mt = sum(mDT[grep('^mt-', GeneName)][[i]])/(sum(mDT[grep('^mt-', GeneName)][[i]]) + sum(hDT[grep('^MT-', GeneName)][[i]])) *100
# print(paste0(i,' Spike-in % by nuclear genes: ',si_nuc))
# print(paste0(i,' Spike-in % by mito genes: ',si_mt))}
# 
# for (i in x){
# sum_all = sum(mDT[[i]])
# sum_mt = sum(mDT[grep('^mt-', GeneName)][[i]])
# print(paste0(i,' Spike-in read number all transcripts: ',sum_all))
# print(paste0(i,' Spike-in read number mito transcripts: ',sum_mt))
# print(paste0(i,' Spike-in read number mito transcripts inc ncRNA: ',sum_mt))}

# 
# 
# ######## Calculate RPK, RPKM, RPKS and make tables ######### 


pernumber=1000
DT_RPK <- copy(mergeDT)
RPKcols <- list()
for (i in c(1:SampNum)) {
	RPKcols <- list.append(RPKcols, round(mergeDT[[samplenames[i]]]/mergeDT$Length*pernumber, digits=2))
	}
DT_RPK[, (samplenames) := RPKcols]

pernumber=1000000	
DT_RPKM <- copy(DT_RPK)
RPKMcols <- list()
for (i in c(1:SampNum)) {
	RPKMcols <- list.append(RPKMcols, round(DT_RPK[[samplenames[i]]]/sum(mergeDT[[samplenames[i]]])*pernumber, digits=2))
	}
DT_RPKM[, (samplenames) := RPKMcols]

pernumber=1000000	
DT_tpm <- copy(DT_RPK)
tpmcols <- list()
for (i in c(1:SampNum)) {
	tpmcols <- list.append(tpmcols, round(DT_RPK[[samplenames[i]]]/sum(DT_RPK[[samplenames[i]]])*pernumber, digits=2))
	}
DT_tpm[, (samplenames) := tpmcols]


# 
# 
# #Write files
# 
if (pseudo == 'include' & ncRNA == 'include') {

write.table(mergeDT, file=paste0(path, fileName, '_reads.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPK, file=paste0(path, fileName, '_RPK.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPKM, file=paste0(path, fileName, '_RPKM.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_tpm, file=paste0(path, fileName, '_tpm.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

}


if (pseudo == 'exclude' & ncRNA == 'include') {

write.table(mergeDT, file=paste0(path, fileName, '_reads_noPseudo.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPK, file=paste0(path, fileName, '_RPK_noPseudo.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPKM, file=paste0(path, fileName, '_RPKM_noPseudo.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_tpm, file=paste0(path, fileName, '_tpm_noPseudo.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

}

if (pseudo == 'exclude' & ncRNA == 'exclude') {

write.table(mergeDT, file=paste0(path, fileName, '_reads_noPseudo_ncRNA.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPK, file=paste0(path, fileName, '_RPK_noPseudo_ncRNA.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_RPKM, file=paste0(path, fileName, '_RPKM_noPseudo_ncRNA.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

write.table(DT_tpm, file=paste0(path, fileName, '_tpm_noPseudo_ncRNA.txt'), row.names=FALSE, col.names=c('GeneName',samplenames,'Length', 'GeneID'), sep=("\t"), quote=FALSE)

}

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/AddGeneName_RPKMS_customGTF.R')
