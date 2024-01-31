library(data.table)
library('scales')
library('vioplot')

wid = 3.5 # width of plot 2.4
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2

sets=c("WT","LRP")

# Set up vectors to hold normalized relative production vals
nRNR2_1s=list()
nMito_1s=list()
nNuc_1s=list()
nRNR2_2s=list()
nMito_2s=list()
nNuc_2s=list()


for (set in sets) {
# Abundance
print('Abundance')
type = 'RPKM' # reads, RPK, RPKM, RPKS 
HLcorrected = 'no'

if (HLcorrected == 'yes') {
time = '1592'
addon = 'HLcorr'
}
if (HLcorrected == 'no') {
time = '100000'
addon = 'notHLcorr'
}


folderandmethod = 'TL10_HEK_2023_11/RPK/TL10_t5MMinformed5' # 
RNA1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
folderandmethod = 'HEK_TL10_LRPPRC_2022_07/RPK/TL10_t5MTMMinformed6' # 
RNA1mito <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))



folderandmethod = 'TL12_HEK_2023_11/RPK/TL12_t5MMinformed5' # Hela_TL4_combined/RPK/TL4_t5MTMMinformed6
RNA2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
folderandmethod = 'HEK_TL12_LRPPRC_2022_08/RPK/TL12_t5MTMMinformed6' # 
RNA2mito <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/', folderandmethod,'_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))






mitogenes = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6') #, 'MT-ND6', '7S')

nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

nucHKgenes = c('GAPDH',  'C1orf43',  'CHMP2A', 'EMC7', 'GPI', 'PSMB2', 'PSMB4',  'MT-CO2',  'RAB7A', 'REEP5', 'SNRPD3','VCP','VPS29') # except for GAPDH, source here: https://www.cell.com/trends/genetics/pdf/S0168-9525%2813%2900089-9.pdf

Nuc_1=NULL
Nuc_2=NULL

if (set == 'WT') {c=2}
if (set == 'LRP') {c=3}

RNR2_1 = RNA1mito[GeneName == 'MT-RNR2'][,..c]
Mito_1 = RNA1mito[GeneName %in% mitogenes][,..c]
Nuc_1 = RNA1[GeneName %in% nucgenes][,..c]
NucHK_1 = RNA1[GeneName %in% nucHKgenes][,..c]
RNR2_2 = RNA2mito[GeneName == 'MT-RNR2'][,..c]
Mito_2 = RNA2mito[GeneName %in% mitogenes][,..c]
Nuc_2 = RNA2[GeneName %in% nucgenes][,..c]
NucHK_2 = RNA2[GeneName %in% nucHKgenes][,..c]

# Remove 0 values
Nuc_1 <- Nuc_1[Nuc_1[[1]]!=0]
Nuc_2 <- Nuc_2[Nuc_2[[1]]!=0]
NucHK_1 <- NucHK_1[NucHK_1[[1]]!=0]
NucHK_2 <- NucHK_2[NucHK_2[[1]]!=0]


# write.table(Mito_1, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(Mito_2, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(Mito_3, '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_mito_Ab_2023_09.txt', row.names=F, quote=F)
# write.table(setorder(RNA1[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(RNA2[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(RNA3[GeneName %in% nucgenes][,1:2], cols='GeneName'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_nuc_Ab_2023_09.txt', row.names=F, quote=F, sep='\t')

catnames = c('RNR2_1','Mito_1','NucOX_1','NucHK_1','RNR_2','Mito_2','NucOX_2','NucHK_2')

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/AbundanceBoxplot_HEK_',set,'_',type,'_2023_011.pdf'), width = wid, height = ht)

vec1=c(RNR2_1,Mito_1, Nuc_1, NucHK_1,RNR2_2,Mito_2, Nuc_2, NucHK_2)

xx=boxplot(vec1, ylab = paste0('Abundance (',type, ')'), main = paste0(set,' Abundance mRNA'), xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)

# Add in points
# points(jitter(xx$group), xx$out, pch=16, cex=.5, col='black') #outliers only
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(NucHK_1)), rep(5,nrow(RNR2_2)), rep(6,nrow(Mito_2)), rep(7,nrow(Nuc_2)), rep(8,nrow(NucHK_2))) 

points(jitter(pointsvec), unlist(vec1), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/Abundance_HEK_',set,'_rRNA_mito_nuc.txt')))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()


# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/AbundanceBoxplot_HEK_',set,'_rep1_',type,'_2023_11.pdf'), width = wid*.66, height = ht)

boxplot(vec1[1:3], ylab = paste0('Abundance (',type, ')'), main = paste0(set,' Abundance mRNA'), xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()





# Half Life
print('Half life')

if (set == 'WT') {method = 'PcMTnorRNA'}
if (set == 'LRP') {method = 'PcBrute1'}


HLnuc1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/TL10_HEK_2023_11/HalfLife/TL10_',set,'_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


HLmito1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/HalfLife/TL10_',set,'_MT_t5MTMMinformed6_modeAll_',method,'_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


HLnuc2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/TL12_HEK_2023_11/HalfLife/TL12_',set,'_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

HLmito2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/HalfLife/TL12_',set,'_MT_t5MTMMinformed6_modeAll_',method,'_FracNew_halflives_corr_',time,'min_v2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


RNR2_1 = HLmito1[Gene == 'MT-RNR2'][,3]
Mito_1 = HLmito1[Gene %in% mitogenes][,3]
Nuc_1 = HLnuc1[Gene %in% nucgenes][,3]
NucHK_1 = HLnuc1[Gene %in% nucHKgenes][,3]
RNR2_2 = HLmito2[Gene == 'MT-RNR2'][,3]
Mito_2 = HLmito2[Gene %in% mitogenes][,3]
Nuc_2 = HLnuc2[Gene %in% nucgenes][,3]
NucHK_2 = HLnuc2[Gene %in% nucHKgenes][,3]

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/HLBoxplot_HEK_',set,'_2023_11_',addon,'.pdf'), width = wid, height = ht)

vec2=c(RNR2_1,Mito_1, Nuc_1,NucHK_1,RNR2_2,Mito_2, Nuc_2, NucHK_2)

xx=boxplot(vec2, ylab = 'Half-life (min)', main = paste0(set,' Turnover mRNA'), xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)

# Add in points
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(NucHK_1)), rep(5,nrow(RNR2_2)), rep(6,nrow(Mito_2)), rep(7,nrow(Nuc_2)), rep(8,nrow(NucHK_2))) 

points(jitter(pointsvec), unlist(vec2), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/Halflife_HEK_',set,'_rRNA_mito_nuc_',addon,'.txt')))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()

# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/HLBoxplot_HEK_',set,'_rep1_2023_11_',addon,'.pdf'), width = wid*.66, height = ht)
boxplot(vec2[1:3], ylab = 'Half-life (min)', main = paste0(set,' Turnover mRNA'), xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()

# write.table(setorder(HLnuc1[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL3_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(HLnuc2[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL4_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')
# write.table(setorder(HLnuc3[Gene %in% nucgenes], cols='Gene'), '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL1_nuc_HL_2023_09.txt', row.names=F, quote=F, sep='\t')




# Relative production mRNA
print('Relative production')

HL1 = rbind(HLnuc1, HLmito1)
HL2 = rbind(HLnuc2, HLmito2)


RNA_HL1nucHK = merge(RNA1[GeneName %in% nucHKgenes][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1nuc = merge(RNA1[GeneName %in% nucgenes][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1rrna = merge(RNA1mito[GeneName == 'MT-RNR2'][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1mito = merge(RNA1mito[GeneName %in% mitogenes][,c(1:2)], HL1[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL1 = rbind(RNA_HL1nucHK,RNA_HL1nuc, RNA_HL1rrna, RNA_HL1mito)

RNA_HL2nucHK = merge(RNA2[GeneName %in% nucHKgenes][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2nuc = merge(RNA2[GeneName %in% nucgenes][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2rrna = merge(RNA2mito[GeneName == 'MT-RNR2'][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2mito = merge(RNA2mito[GeneName %in% mitogenes][,c(1:2)], HL2[, c('Gene', 'Half.Life')], by.x='GeneName', by.y='Gene')
RNA_HL2 = rbind(RNA_HL2nucHK,RNA_HL2nuc, RNA_HL2rrna, RNA_HL2mito)



setnames(RNA_HL1, 2, 'Abundance')
setnames(RNA_HL2, 2, 'Abundance')


RNA_HL1[, RelProd := Abundance/Half.Life]
RNA_HL2[, RelProd := Abundance/Half.Life]


RNR2_1 = RNA_HL1[GeneName == 'MT-RNR2'][,4]
Mito_1 = RNA_HL1[GeneName %in% mitogenes][,4]
Nuc_1 = RNA_HL1[GeneName %in% nucgenes][,4]
NucHK_1 = RNA_HL1[GeneName %in% nucHKgenes][,4]
RNR2_2 = RNA_HL2[GeneName == 'MT-RNR2'][,4]
Mito_2 = RNA_HL2[GeneName %in% mitogenes][,4]
Nuc_2 = RNA_HL2[GeneName %in% nucgenes][,4]
NucHK_2 = RNA_HL2[GeneName %in% nucHKgenes][,4]

# Normalize to NucHK

nRNR2_1=RNR2_1[[1]]/median(NucHK_1[[1]])
nMito_1=Mito_1[[1]]/median(NucHK_1[[1]])
nNuc_1=Nuc_1[[1]]/median(NucHK_1[[1]])
nRNR2_2=RNR2_2[[1]]/median(NucHK_2[[1]])
nMito_2=Mito_2[[1]]/median(NucHK_2[[1]])
nNuc_2=Nuc_2[[1]]/median(NucHK_2[[1]])

nRNR2_1s[[length(nRNR2_1s)+1]] = nRNR2_1
nMito_1s[[length(nMito_1s)+1]] = nMito_1
nNuc_1s[[length(nNuc_1s)+1]] = nNuc_1
nRNR2_2s[[length(nRNR2_2s)+1]] = nRNR2_2
nMito_2s[[length(nMito_2s)+1]] = nMito_2
nNuc_2s[[length(nNuc_2s)+1]] = nNuc_2

# Store values for plotting WT and LRP together




pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot_HEK_',set,'_2023_11_',addon,'.pdf'),width = wid, height = ht)

vec3=c(RNR2_1,Mito_1, Nuc_1, NucHK_1,RNR2_2,Mito_2, Nuc_2, NucHK_2)

xx=boxplot(vec3, ylab = 'Production rate', main = paste0(set,' Relative production mRNA'), xlab = '', names=catnames, pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
# Add in points
pointsvec=c(rep(1,nrow(RNR2_1)), rep(2,nrow(Mito_1)), rep(3,nrow(Nuc_1)), rep(4,nrow(NucHK_1)), rep(5,nrow(RNR2_2)), rep(6,nrow(Mito_2)), rep(7,nrow(Nuc_2)), rep(7,nrow(NucHK_2))) 

points(jitter(pointsvec), unlist(vec3), pch=16, cex=.4, col=alpha('dodgerblue', .5))

dev.off()

sink(file.path(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/BoxplotData/RelProd_HEK_',set,'_rRNA_mito_nuc_',addon,'.txt')))
print('Boxplot output', quote=FALSE)
print(catnames, quote=FALSE)
print(xx)
sink()

# And plot just rep1
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot_HEK_',set,'_rep1_2023_11_',addon,'.pdf'),width = wid*.66, height = ht)
boxplot(vec3[1:3], ylab = 'Production rate', main = paste0(set,' Relative production mRNA'), xlab = '', names=catnames[1:3], pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)
dev.off()


write.table(setorder(RNA_HL1, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL10_HEK_',set,'_ab_HL_relprod_2023_11_',addon,'.txt'), row.names=F, quote=F, sep='\t')
write.table(setorder(RNA_HL2, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/TL12_HEK_',set,'_ab_HL_relprod_2023_11_',addon,'.txt'), row.names=F, quote=F, sep='\t')

}


# Plot WT and LRP boxplots (normalized to NucHK) side by side


vec4=c(nRNR2_1s, nMito_1s,nNuc_1s, nRNR2_2s, nMito_2s, nNuc_2s)

catnames2 = rep(c('RNR2_1','Mito_1','NucOX_1','RNR_2','Mito_2','NucOX_2'), each=2)

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/RelProdBoxplot_HEK_both_2023_11_',addon,'.pdf'),width = wid, height = ht*1.25)

xx=boxplot(vec4, ylab = "Production rate norm to HK genes' med", main = "Relative production mRNA", xlab = "", names=catnames2, pch=16, outcex=.2, col=c("purple", "orange"), outcol=alpha("black", .2), cex.axis=txt, log="y",outline=FALSE, las=las) #,legend.text=c("WT", "LRPPRC KO"), args.legend = list(x='bottomleft', bty = 'n', fill = c("purple", "orange"), cex = .8))

legend("bottomleft", c("WT", "LRPPRC KO"), bty = 'n', fill = c("purple", "orange"), cex = .6)

dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots_LRPPRC.r')

