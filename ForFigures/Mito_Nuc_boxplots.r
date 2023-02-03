library(data.table)
library('scales')


wid = 2.4 # width of plot
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2
# Abundance

type = 'RPKM' # reads, RPK, RPKM, RPKS 

RNA1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/RPK/TL3_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

RNA2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/RPK/TL4_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


mitogenes = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6')

nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

RNR2_1 = RNA1[GeneName == 'MT-RNR2'][,2]
Mito_1 = RNA1[GeneName %in% mitogenes][,2]
Nuc_1 = RNA1[GeneName %in% nucgenes][,2]
# Remove 0 values
Nuc_1 <- Nuc_1[Nuc_1[[1]]!=0]


catnames = c('RNR2','Mito', 'Nuc')

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/AbundanceBoxplot_',type,'.pdf'), width = wid, height = ht)

xx=boxplot(c(RNR2_1, Mito_1, Nuc_1), ylab = paste0('Abundance (',type, ')'), main = 'Abundance mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)


dev.off()


# Half Life


HLnuc1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_Nuc_t5MMinformed5_lenient_modeAll_Mitogenes_FracNew_halflives_corr_1592min.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_1592min.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLnuc2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_Nuc_t5MMinformed5_lenient_modeAll_Mitogenes_FracNew_halflives_corr_1592min.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_1592min.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

RNR2_1 = HLmito1[Gene == 'MT-RNR2']$Half.Life
Mito_1 = HLmito1[Gene %in% mitogenes]$Half.Life
Nuc_1 = HLnuc1$Half.Life

catnames = c('RNR2','Mito', 'Nuc')

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/HLBoxplot.pdf'), width = wid, height = ht)

xx=boxplot(RNR2_1, Mito_1, Nuc_1, ylab = 'Half-life (min)', main = 'Turnover mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)


dev.off()



# Relative production mRNA

HL1 = rbind(HLnuc1, HLmito1)
HL2 = rbind(HLnuc2, HLmito2)


RNA_HL1 = merge(RNA1[,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')

setnames(RNA_HL1, 2, 'Abundance')

RNA_HL1[, RelProd := Abundance/Half.Life]

RNR2_1 = RNA_HL1[GeneName == 'MT-RNR2']$RelProd
Mito_1 = RNA_HL1[GeneName %in% mitogenes]$RelProd
Nuc_1 = RNA_HL1[!grep('MT-', GeneName)]$RelProd


pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot.pdf'),width = wid, height = ht)

xx=boxplot(RNR2_1, Mito_1, Nuc_1, ylab = 'Production rate', main = 'Relative production mRNA', xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)


dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots.r')

