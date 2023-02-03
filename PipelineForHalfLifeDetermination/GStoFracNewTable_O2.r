library(data.table)
library('scales')
# library('calibrate')

### USE (run from same directory that GS directory is in
### sbatch -p short -t 0-01:00 --wrap="Rscript /n/groups/churchman/mc348/TimelapseSeq/Scripts/GStoFracNewTable_O2.r"

args <- commandArgs(trailingOnly = TRUE)

measure = 'MAP' # MAP mean
#tps = c('0', '15', '30','45', '60', '90', '120', '240')
# tps = c('0', '7', '15', '30', '45','60', '90', '120', '240')
# tps = c('0', '15','30', '60')
#tps = c('0', '30', '60', '120')
#tps = c('0', '60', '120', '240')

tps = c('0', '30', '60')
#suffixes = c('m', 'm_4sU', 'm_4sU6sG_200')
#suffixes = c('m', 'm_4sU', 'm_4sU6sG')
suffixes = c('m')

# args <- commandArgs(trailingOnly = TRUE)
Expt = args[1] #'TL13'
expt = args[2]  # longer version, in file name
Set = args[3] #'MTdsrRNA' # Nuc MT MTdsrRNA
Modifier = args[4] #'t5MTMMinformed6_modeAll' # t5MMinformed5_lenient_modeAll
# TL1_MT_t5MMinformed4_modeAll TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA t5MMinformed5_lenient_modeAll
Pcset=args[5]

path = paste0(getwd(), '/')

if (Set != 'Nuc') {GSpath = paste0('GS_',expt,'_',Set,'_',Modifier, Pcset)}
if (Set == 'Nuc') {GSpath = paste0('GS_',expt,'_',Set,'_',Modifier)}


filename = paste0(expt,'_',Set, '_', Modifier) 

DT = data.table(read.table(paste0(path,GSpath,'/',filename,'.tsv'), sep='\t', header=TRUE))
# Readcount filter
# rcCOI <- paste0(expt, '_', tps, 'm.Readcount')

#COI <- paste0(Expt, '_', tps, 'm.', measure) 
COI <- paste0(Expt, '_', tps,suffixes, '.', measure)
DT <- cbind(DT[, list(Symbol)], DT[,..COI])
DT <- DT[!grep('MT-T', DT$Symbol)] 

# # matplot(t(na.omit(DT[DT$NHC2_240m.MAP != 0 & DT$NHC2_120m.MAP != 0 & DT$NHC2_90m.MAP != 0  & DT$NHC2_60m.MAP != 0])), type='l')


if (Set != 'Nuc') {write.table(DT, file=paste0(path,'FracNew/',filename,Pcset, '_FracNew.txt'), row.names=FALSE, col.names = c('Symbol',tps),sep=("\t"), quote=FALSE)}
if (Set == 'Nuc') {write.table(DT, file=paste0(path,'FracNew/',filename, '_FracNew.txt'), row.names=FALSE, col.names = c('Symbol',tps),sep=("\t"), quote=FALSE)}

# If it's the nuclear set also filter for nuclear mito genes only

if (Set == 'Nuc') {
GOI=data.table(read.table(paste0(path, '../SeqFiles/Mito genes_newSymbols.txt'), sep='\t', header=TRUE))

# Filter DT
filtDT <- DT[DT$Symbol %in% GOI$Mito.genes]

write.table(filtDT, file=paste0(path,'FracNew/',filename, '_Mitogenes_FracNew.txt'), row.names=FALSE, col.names = c('Symbol',tps),sep=("\t"), quote=FALSE)

}






