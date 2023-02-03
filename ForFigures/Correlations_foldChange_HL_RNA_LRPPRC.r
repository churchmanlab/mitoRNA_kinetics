library('scales')
library(data.table)
library(rlist)
library(Rfast)

HLwt1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/HalfLife/TL10_WT_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt1', 'HLnorm')))

HLmut1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/HalfLife/TL10_LRP_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_30minEdit_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut1', 'HLnorm')))

HLwt2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/HalfLife/TL12_WT_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt2', 'HLnorm')))

HLmut2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/HalfLife/TL12_LRP_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut2', 'HLnorm')))

# Merge because genes are not in the same order above
HLs_tmp1 <- merge(HLwt1[, c('Gene', 'HL_wt1')], HLmut1[, c('Gene', 'HL_mut1')], by='Gene')
HLs_tmp2 <- merge(HLs_tmp1[, c('Gene', 'HL_wt1', 'HL_mut1')], HLwt2[, c('Gene', 'HL_wt2')], by='Gene')
HLs <- merge(HLs_tmp2[, c('Gene', 'HL_wt1', 'HL_mut1', 'HL_wt2')], HLmut2[, c('Gene', 'HL_mut2')], by='Gene')
# Get fold change
HLs[, HL_FC1 := log(HL_mut1/HL_wt1, 2)]
HLs[, HL_FC2 := log(HL_mut2/HL_wt2, 2)]

# Remove RNR1
HLs <- HLs[Gene != 'MT-RNR1']


# RNA abundance
RNA1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/RPK/TL10_t5MTMMinformed6_featureCounts_multi_RPK_noPseudo.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))

RNA2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/RPK/TL12_t5MTMMinformed6_featureCounts_multi_RPK_noPseudo.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))

# Use RNA read count table without ncRNA for normalizing
RNAfornorm1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/RPK/TL10_t5MTMMinformed6_featureCounts_multi_reads_noPseudo_ncRNA.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))

RNAfornorm2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/RPK/TL12_t5MTMMinformed6_featureCounts_multi_reads_noPseudo_ncRNA.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))

# Normalize to counts on nuclear genes
# Remove RNR1 and ND6, tRNAs, anti genes
RNAfornorm1 <- RNAfornorm1[GeneName != '7S']
RNAfornorm1 <- RNAfornorm1[!grep('^MT-', GeneName)]
RNAfornorm2 <- RNAfornorm2[GeneName != '7S']
RNAfornorm2 <- RNAfornorm2[!grep('^MT-', GeneName)]

# Rep1
SampNum=ncol(RNA1)-3
samplenames = colnames(RNA1)[2:(SampNum+1)]

pernumber = 1000000
RNA_RPKMnuc1 <- copy(RNA1)
RPKMnuccols <- list()
for (i in c(1:SampNum)) {
	RPKMnuccols <- list.append(RPKMnuccols, round(RNA1[[samplenames[i]]]/sum(RNAfornorm1[[samplenames[i]]])*pernumber, digits=2))
	}
RNA_RPKMnuc1[, (samplenames) := RPKMnuccols]

# Rep2
SampNum=ncol(RNA2)-3
samplenames = colnames(RNA2)[2:(SampNum+1)]

pernumber = 1000000
RNA_RPKMnuc2 <- copy(RNA2)
RPKMnuccols <- list()
for (i in c(1:SampNum)) {
	RPKMnuccols <- list.append(RPKMnuccols, round(RNA2[[samplenames[i]]]/sum(RNAfornorm2[[samplenames[i]]])*pernumber, digits=2))
	}
RNA_RPKMnuc2[, (samplenames) := RPKMnuccols]


# Merge HL and normalized RNA tables
DT_tmp <- merge(HLs, RNA_RPKMnuc1[, c('GeneName','TL10_0m_WT', 'TL10_0m_LRP')], by.x='Gene', by.y='GeneName') 
DT <- merge(DT_tmp, RNA_RPKMnuc2[, c('GeneName','TL12_0m_WT', 'TL12_0m_LRP')], by.x='Gene', by.y='GeneName') 
# Get abundance fold change
DT[, RNA_FC1 := log(TL10_0m_LRP/TL10_0m_WT, 2)]
DT[, RNA_FC2 := log(TL12_0m_LRP/TL12_0m_WT, 2)]


heavyStrand = c('MT-RNR2', 'MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB')
lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1', '7S')


# Heavy strand

DTsub <- DT[DT$Gene %in% heavyStrand]


# First just plot mut vs wt
xsamps = c('TL10_0m_WT', 'TL12_0m_WT')
ysamps = c('TL10_0m_LRP', 'TL12_0m_LRP')

# Get x and y
x = rowMeans(DTsub[,..xsamps])
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# remove RNR2
x = x[1:10]
y = y[1:10]
genes = genes[1:10]

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..ysamps], 1, sd)    
xsd = xsd[1:10]
ysd = ysd[1:10]

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/Abundance_LRPvsWT_heavy.pdf'), width = 3, height = 3.5)

limits=c(-4.1, 2.2)

xName = 'WT abundance'
yName = 'LRPPRC KO abundance'

colptsize = 1

# LRP abundance vs WT
limits=c(0,10000)
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=limits, xlim=limits)
 # 

# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)


# Diagonal line
abline(0,1, lwd = 0.5, col='black', lty=2)
# Trendline
# Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
points(x[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)], y[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)],pch = 21, cex = colptsize, bg='dodgerblue',col = 'black')
points(x[grep('MT-CYB', genes)], y[grep('MT-CYB', genes)],pch = 21, cex = colptsize, bg = 'aquamarine2',col = 'black')
points(x[grep('MT-CO1|MT-CO2|MT-CO3', genes)], y[grep('MT-CO1|MT-CO2|MT-CO3', genes)],pch = 21, cex = colptsize, bg = 'pink3',col = 'black')
points(x[grep('MT-ATP8-6', genes)], y[grep('MT-ATP8-6', genes)],pch = 21, cex = colptsize, bg = 'grey60',col = 'black')

# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# 7S RNA
# points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 16, cex = colptsize, col = 'brown')

text(min(limits)+(max(limits)-min(limits))/4, max(limits)*(4/5), labels = mylabel, cex = .8)

legend("topright", legend=c("CI", "CIII", "CIV", "CV"), cex=.8, bty="n", text.col =  c('dodgerblue', 'aquamarine2', 'pink3','grey60'),  ncol=1) 

dev.off()









DTsub <- DT[DT$Gene %in% heavyStrand]

# Now mut/wt FC vs wt HL
xsamps = c('HL_wt1', 'HL_wt2')
ysamps = c('RNA_FC1', 'RNA_FC2')

# Get x and y
x = rowMeans(DTsub[,..xsamps])
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# remove RNR2
x = x[1:10]
y = y[1:10]
genes = genes[1:10]

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..ysamps], 1, sd)    
xsd = xsd[1:10]
ysd = ysd[1:10]

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/LRPeffect_vs_HL_heavy.pdf'), width = 3, height = 3.5)


xName = 'WT half-life'
yName = 'log2 LRPPRC KO / WT'

colptsize = 1

# LRP abundance vs WT
ylimits=c(-4.1, -.25)
xlimits=c(10, 140)
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=ylimits, xlim=xlimits)
 # , ylim=limits, xlim=limits)

# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)


# Diagonal line
abline(0,1, lwd = 0.5, col='black', lty=2)
# Trendline
# Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
points(x[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)], y[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)],pch = 21, cex = colptsize, bg='dodgerblue',col = 'black')
points(x[grep('MT-CYB', genes)], y[grep('MT-CYB', genes)],pch = 21, cex = colptsize, bg = 'aquamarine2',col = 'black')
points(x[grep('MT-CO1|MT-CO2|MT-CO3', genes)], y[grep('MT-CO1|MT-CO2|MT-CO3', genes)],pch = 21, cex = colptsize, bg = 'pink3',col = 'black')
points(x[grep('MT-ATP8-6', genes)], y[grep('MT-ATP8-6', genes)],pch = 21, cex = colptsize, bg = 'grey60',col = 'black')

# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# 7S RNA
# points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 16, cex = colptsize, col = 'brown')

text(min(xlimits)+(max(xlimits)-min(xlimits))/4, -.5, labels = mylabel, cex = .8)

legend("topright", legend=c("CI", "CIII", "CIV", "CV"), cex=.8, bty="n", text.col =  c('dodgerblue', 'aquamarine2', 'pink3','grey60'),  ncol=1) 

dev.off()









# Now plot RNA fc vs HL fc
DTsub <- DT[DT$Gene %in% heavyStrand]

xsamps = c('HL_FC1', 'HL_FC2')
ysamps = c('RNA_FC1', 'RNA_FC2')

# Get x and y
x = rowMeans(DTsub[,..xsamps]) # HL fold-change
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..xsamps], 1, sd)    

# Get ranges
xrange = rowMaxs(as.matrix(DTsub[,..xsamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..xsamps]), value=TRUE)

yrange = rowMaxs(as.matrix(DTsub[,..ysamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..ysamps]), value=TRUE)



# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/AbundanceFC_vs_HL_FC_heavy.pdf'), width = 3, height = 3.5)

limits=c(-4.1, 2.2)

xName = 'Half-life Log2 FC'
yName = 'Abundance Log2 FC'

colptsize = 1

# Half lives vs abundance
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=limits, xlim=limits)


# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)

# Range bars
# arrows(x-xrange/2, y, x+(xrange/2), y, length=0.05, angle=90, code=3, lwd=.2)
# arrows(x, y-yrange/2, x, y+(yrange/2), length=0.05, angle=90, code=3, lwd=.2)

# Diagonal line
abline(0,1, lwd = 0.5, col='black', lty=2)
abline(h=0, lwd = 0.5)
abline(v=0, lwd = 0.5)
# Trendline
# Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
points(x[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)], y[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)],pch = 21, cex = colptsize, bg='dodgerblue',col = 'black')
points(x[grep('MT-CYB', genes)], y[grep('MT-CYB', genes)],pch = 21, cex = colptsize, bg = 'aquamarine2',col = 'black')
points(x[grep('MT-CO1|MT-CO2|MT-CO3', genes)], y[grep('MT-CO1|MT-CO2|MT-CO3', genes)],pch = 21, cex = colptsize, bg = 'pink3',col = 'black')
points(x[grep('MT-ATP8-6', genes)], y[grep('MT-ATP8-6', genes)],pch = 21, cex = colptsize, bg = 'grey60',col = 'black')

# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# 7S RNA
# points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 16, cex = colptsize, col = 'brown')

text(min(limits)+(max(limits)-min(limits))/4, max(limits)*(4/5), labels = mylabel, cex = .8)

legend("bottomright", legend=c("RNR2","CI", "CIII", "CIV", "CV"), cex=.8, bty="n", text.col =  c('purple','dodgerblue', 'aquamarine2', 'pink3','grey60'),  ncol=1) 

dev.off()



# And a barplot with fold-change

AbFC=y
# Order by most to least
ord=order(AbFC)
AbFC=AbFC[ord]
genesord=genes[ord]
AbFCsdv = ysd[ord]


# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/AbundanceFC_barplot.pdf'), width = 3.5, height = 3.5)

xx=barplot(AbFC, ylim=c(-4.5,1), names=genesord, las=2, cex.names=.7, ylab='LRPPRC KO/WT RPKM (log2)')

arrows(xx, AbFC-AbFCsdv/2, xx, AbFC+(AbFCsdv/2), length=0.05, angle=90, code=3, lwd=.2)


dev.off()





# Light strand

DTsub <- DT[DT$Gene %in% lightStrand]

xsamps = c('HL_FC1', 'HL_FC2')
ysamps = c('RNA_FC1', 'RNA_FC2')

# Get x and y
x = rowMeans(DTsub[,..xsamps]) # HL fold-change
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..xsamps], 1, sd)    



# Get ranges
xrange = rowMaxs(as.matrix(DTsub[,..xsamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..xsamps]), value=TRUE)

yrange = rowMaxs(as.matrix(DTsub[,..ysamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..ysamps]), value=TRUE)


# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/AbundanceFC_vs_HL_FC_light.pdf'), width = 3, height = 3.5)

# limits=c(-4, 1.1)
# ylimits = c(-2,2)

xName = 'Half-life Log2 FC'
yName = 'Abundance Log2 FC'

colptsize = 1

# Half lives vs abundance
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=limits, xlim=limits)

# Range bars
# arrows(x-xrange/2, y, x+(xrange/2), y, length=0.05, angle=90, code=3, lwd=.2)
# arrows(x, y-yrange/2, x, y+(yrange/2), length=0.05, angle=90, code=3, lwd=.2)
# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)

# Diagonal line
abline(0,1, lwd = 0.5, col='black', lty=2)
abline(h=0, lwd = 0.5)
abline(v=0, lwd = 0.5)
# Trendline
# Mito genes
points(x[genes %in% lightStrand], y[genes %in% lightStrand],pch = 17, cex = colptsize, col = 'red')
# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 17, cex = colptsize, col = 'purple')
# 7S RNA
points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 17, cex = colptsize, col = 'brown')

text(min(limits)+(max(limits)-min(limits))/4, max(limits)*(4/5), labels = mylabel, cex = .8)

legend("bottomright", legend=c("heavy", "light", "MT-RNR2", "7S"), cex=.8, bty="n", text.col =  c('orange', 'red', 'purple','brown'),  ncol=1) 

dev.off()




# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_foldChange_HL_RNA_LRPPRC.R')

