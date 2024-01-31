library(data.table)
library('scales')
library('vioplot')

wid = 4 # width of plot 2.4
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2


# Abundance
print('Abundance')
type = 'RPKM' # reads, RPK, RPKM, RPKS 
HLcorrected = 'yes'

if (HLcorrected == 'yes') {
time = '1592'
addon = 'HLcorr'

}
if (HLcorrected == 'no') {
time = '100000'
addon = 'notHLcorr'
}



folderandmethod = 'TL3_Hela_2023_09/RPK/TL3_t5MMinformed5' # Hela_TL3_2020_09/RPK/TL3_t5MTMMinformed6
RNA1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
folderandmethod = 'Hela_TL3_2020_09/RPK/TL3_t5MTMMinformed6' # 
RNA1mito <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

folderandmethod = 'TL4_Hela_2023_09/RPK/TL4_t5MMinformed5' # Hela_TL4_combined/RPK/TL4_t5MTMMinformed6
RNA2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
folderandmethod = 'Hela_TL4_combined/RPK/TL4_t5MTMMinformed6' # 
RNA2mito <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

folderandmethod = 'TL1_Hela_2023_09/RPK/TL1_t5MMinformed4' # Hela_TL1_2020_02/RPK/TL1_t5MMinformed4
RNA3 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
folderandmethod = 'Hela_TL1_2020_02/RPK/TL1_t5MMinformed4' # 
RNA3mito <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


mitogenes = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6') #, 'MT-ND6', '7S')

nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

Nuc_1=NULL
Nuc_2=NULL
Nuc_3=NULL
RNR2_1 = RNA1mito[GeneName == 'MT-RNR2'][,2]
Mito_1 = RNA1mito[GeneName %in% mitogenes][,2]
Nuc_1 = RNA1[GeneName %in% nucgenes][,2]
RNR2_2 = RNA2mito[GeneName == 'MT-RNR2'][,2]
Mito_2 = RNA2mito[GeneName %in% mitogenes][,2]
Nuc_2 = RNA2[GeneName %in% nucgenes][,2]
RNR2_3 = RNA3mito[GeneName == 'MT-RNR2'][,2]
Mito_3 = RNA3mito[GeneName %in% mitogenes][,2]
Nuc_3 = RNA3[GeneName %in% nucgenes][,2]
# Remove 0 values
Nuc_1 <- Nuc_1[Nuc_1[[1]]!=0]
Nuc_2 <- Nuc_2[Nuc_2[[1]]!=0]
Nuc_3 <- Nuc_3[Nuc_3[[1]]!=0]

# write.table(Mito_1, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(Mito_2, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(Mito_3, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(setorder(RNA1[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(RNA2[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(RNA3[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')

catnames = c('RNR2_low_1','Mito_low_1','Nuc_low_1','RNR_low_2','Mito_low_2','Nuc_low_2', 'RNR2_high', 'Mito_high', 'Nuc_high')

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/AbundanceBoxplot_',type,'_2023_09.pdf'), width = wid, height = ht)

vec1=c(RNR2_1,Mito_1, Nuc_1,RNR2_2,Mito_2, Nuc_2, RNR2_3, Mito_3, Nuc_3)

xx=boxplot(vec1, ylab = paste0('Abundance (',type, ')'), main = 'Abundance mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)

# Add in points
# points(jitter(xx$group), xx$out, pch=16, cex=.5, col='black') #outliers only
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(RNR2_2)), rep(5,nrow(Mito_2)), rep(6,nrow(Nuc_2)), rep(7,nrow(RNR2_3)), rep(8,nrow(Mito_3)), rep(9,nrow(Nuc_3))) 

points(jitter(pointsvec), unlist(vec1), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/Abundance_rRNA_mito_nuc.txt'))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()


# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/AbundanceBoxplot_rep1_',type,'_2023_09.pdf'), width = wid*.66, height = ht)

boxplot(vec1[1:3], ylab = paste0('Abundance (',type, ')'), main = 'Abundance mRNA', xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()





# Half Life
print('Half life')

if (HLcorrected == 'yes') {ver = 'min'}
if (HLcorrected == 'no') {ver = 'min_v2'}

HLnuc1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/TL3_Hela_2023_09/HalfLife/TL3_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_',time,ver,'.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLnuc2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/TL4_Hela_2023_09/HalfLife/TL4_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_',time,ver,'.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLnuc3 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/TL1_Hela_2023_09/HalfLife/TL1_Nuc_t5MMinformed4_lenient_modeAll_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito3 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL1_2020_02/HalfLife/TL1_MT_t5MMinformed4_modeAll_FracNew_halflives_corr_',time,ver,'.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

RNR2_1 = HLmito1[Gene == 'MT-RNR2'][,3]
Mito_1 = HLmito1[Gene %in% mitogenes][,3]
Nuc_1 = HLnuc1[Gene %in% nucgenes][,3]
RNR2_2 = HLmito2[Gene == 'MT-RNR2'][,3]
Mito_2 = HLmito2[Gene %in% mitogenes][,3]
Nuc_2 = HLnuc2[Gene %in% nucgenes][,3]
RNR2_3 = HLmito3[Gene == 'MT-RNR2'][,3]
Mito_3 = HLmito3[Gene %in% mitogenes][,3]
Nuc_3 = HLnuc3[Gene %in% nucgenes][,3]

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/HLBoxplot_2023_09_',addon,'.pdf'), width = wid, height = ht)

vec2=c(RNR2_1,Mito_1, Nuc_1,RNR2_2,Mito_2, Nuc_2, RNR2_3, Mito_3, Nuc_3)

xx=boxplot(vec2, ylab = 'Half-life (min)', main = 'Turnover mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)

# Add in points
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(RNR2_2)), rep(5,nrow(Mito_2)), rep(6,nrow(Nuc_2)), rep(7,nrow(RNR2_3)), rep(8,nrow(Mito_3)), rep(9,nrow(Nuc_3))) 

points(jitter(pointsvec), unlist(vec2), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/Halflife_rRNA_mito_nuc_',addon,'.txt')))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()

# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/HLBoxplot_rep1_2023_09_',addon,'.pdf'), width = wid*.66, height = ht)
boxplot(vec2[1:3], ylab = 'Half-life (min)', main = 'Turnover mRNA', xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()

# write.table(setorder(HLnuc1[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(HLnuc2[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(HLnuc3[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')




# Relative production mRNA
print('Relative production')

HL1 = rbind(HLnuc1, HLmito1)
HL2 = rbind(HLnuc2, HLmito2)
HL3 = rbind(HLnuc3, HLmito3)


RNA_HL1nuc = merge(RNA1[GeneName %in% nucgenes][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1rrna = merge(RNA1mito[GeneName == 'MT-RNR2'][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1mito = merge(RNA1mito[GeneName %in% mitogenes][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1 = rbind(RNA_HL1nuc, RNA_HL1rrna, RNA_HL1mito)

RNA_HL2nuc = merge(RNA2[GeneName %in% nucgenes][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2rrna = merge(RNA2mito[GeneName == 'MT-RNR2'][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2mito = merge(RNA2mito[GeneName %in% mitogenes][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2 = rbind(RNA_HL2nuc, RNA_HL2rrna, RNA_HL2mito)

RNA_HL3nuc = merge(RNA3[GeneName %in% nucgenes][,c(1:2)], HL3[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL3rrna = merge(RNA3mito[GeneName == 'MT-RNR2'][,c(1:2)], HL3[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL3mito = merge(RNA3mito[GeneName %in% mitogenes][,c(1:2)], HL3[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL3 = rbind(RNA_HL3nuc, RNA_HL3rrna, RNA_HL3mito)

setnames(RNA_HL1, 2, 'Abundance')
setnames(RNA_HL2, 2, 'Abundance')
setnames(RNA_HL3, 2, 'Abundance')

RNA_HL1[, RelProd := Abundance/Half.Life]
RNA_HL2[, RelProd := Abundance/Half.Life]
RNA_HL3[, RelProd := Abundance/Half.Life]

RNR2_1 = RNA_HL1[GeneName == 'MT-RNR2'][,4]
Mito_1 = RNA_HL1[GeneName %in% mitogenes][,4]
Nuc_1 = RNA_HL1[GeneName %in% nucgenes][,4]
RNR2_2 = RNA_HL2[GeneName == 'MT-RNR2'][,4]
Mito_2 = RNA_HL2[GeneName %in% mitogenes][,4]
Nuc_2 = RNA_HL2[GeneName %in% nucgenes][,4]
RNR2_3 = RNA_HL3[GeneName == 'MT-RNR2'][,4]
Mito_3 = RNA_HL3[GeneName %in% mitogenes][,4]
Nuc_3 = RNA_HL3[GeneName %in% nucgenes][,4]


pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot_2023_09_',addon,'.pdf'),width = wid, height = ht)

vec3=c(RNR2_1,Mito_1, Nuc_1,RNR2_2,Mito_2, Nuc_2, RNR2_3, Mito_3, Nuc_3)

xx=boxplot(vec3, ylab = 'Production rate', main = 'Relative production mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
# Add in points
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(RNR2_2)), rep(5,nrow(Mito_2)), rep(6,nrow(Nuc_2)), rep(7,nrow(RNR2_3)), rep(8,nrow(Mito_3)), rep(9,nrow(Nuc_3))) 

points(jitter(pointsvec), unlist(vec3), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/RelProd_rRNA_mito_nuc_',addon,'.txt')))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()

# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot_rep1_2023_09_',addon,'.pdf'),width = wid*.66, height = ht)
boxplot(vec3[1:3], ylab = 'Production rate', main = 'Relative production mRNA', xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()


write.table(setorder(RNA_HL1, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_ab_HL_relprod_2023_09_',addon,'.txt'), row.names=F, quote=F, sep='\t')
write.table(setorder(RNA_HL2, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_ab_HL_relprod_2023_09_',addon,'.txt'), row.names=F, quote=F, sep='\t')
write.table(setorder(RNA_HL3, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_ab_HL_relprod_2023_09_',addon,'.txt'), row.names=F, quote=F, sep='\t')

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots.r')

