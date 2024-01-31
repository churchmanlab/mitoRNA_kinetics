library('scales')
library(data.table)
library(rlist)
library(Rfast)

directory1 <- 'HEK_TL10_LRPPRC_2022_07' # HEK_TL14_LRPPRC_2023_05
directory2 <- 'HEK_TL12_LRPPRC_2022_08'
MapMethod <- 't5MTMMinformed6'
Expt <- 'TL12'

# WT1 <- 'TL14_WT_A' # TL10_WT
# WT2 <- 'TL14_WT_B' # TL12_WT
# mut1 <- 'TL14_LRP_A'
# mut2 <- 'TL14_LRP_A'

WT1 <- 'TL10_WT' # TL10_WT
WT2 <- 'TL12_WT' # TL12_WT
mut1 <- 'TL10_LRP'
mut2 <- 'TL12_LRP'

WTmethod= '_PcMTnorRNA_' # _PcMTnorRNA_ (This for TL10 and TL12) # _ (TL14)
WTmod='halflives' #  halflives (This for TL10 and TL12) _30minEdit_halflives_ 30minEdit2_halflives # halflives (TL14)

mutmethod='_PcBrute1_' # PcMTnorRNA PcBrute1 PcBrute2 # PcBrute1 uses the RNR2 fraction new calculated from the average half-lives. PcBrute2 directly uses the average RNR2 fractions new # PcBrute1 (TL10, TL12) # _ (TL14)
mutmod='halflives' # halflives (for TL10, TL12, TL14) 30minEdit_halflives 30minEdit2_halflives

filename <- paste0(WT1,'_MT_',MapMethod,'_modeAll',WTmethod,'FracNew_',WTmod,'_corr_100000min_v2.txt')
# TL14_WT_B_MT_t5MTMMinformed6_modeAll_FracNew_halflives_corr_10000min_v2.txt
HLwt1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt1', 'HLnorm')))

filename <- paste0(mut1,'_MT_',MapMethod,'_modeAll',mutmethod,'FracNew_',mutmod,'_corr_100000min_v2.txt')
HLmut1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut1', 'HLnorm')))

filename <- paste0(WT2,'_MT_',MapMethod,'_modeAll',WTmethod,'FracNew_',WTmod,'_corr_100000min_v2.txt')
# TL14_WT_B_MT_t5MTMMinformed6_modeAll_FracNew_halflives_corr_10000min_v2.txt
HLwt2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt2', 'HLnorm')))

filename <- paste0(mut2,'_MT_',MapMethod,'_modeAll',mutmethod,'FracNew_',mutmod,'_corr_100000min_v2.txt')
HLmut2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut2', 'HLnorm')))





# Merge because genes are not in the same order above
HLs_tmp1 <- merge(HLwt1[, c('Gene', 'HL_wt1')], HLmut1[, c('Gene', 'HL_mut1')], by='Gene')
HLs_tmp2 <- merge(HLs_tmp1[, c('Gene', 'HL_wt1', 'HL_mut1')], HLwt2[, c('Gene', 'HL_wt2')], by='Gene')
HLs <- merge(HLs_tmp2[, c('Gene', 'HL_wt1', 'HL_mut1', 'HL_wt2')], HLmut2[, c('Gene', 'HL_mut2')], by='Gene')
# Get fold change
HLs[, HL_FC1 := log(HL_mut1/HL_wt1, 2)]
HLs[, HL_FC2 := log(HL_mut2/HL_wt2, 2)]
HLs[, HL_recFC1 := log(HL_wt1/HL_mut1, 2)]
HLs[, HL_recFC2 := log(HL_wt2/HL_mut2, 2)]

# Remove RNR1
HLs <- HLs[Gene != 'MT-RNR1']


# RNA abundance
if (Expt == 'TL14') {
RNA <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/RPK/',Expt,'_',MapMethod,'_featureCounts_multi_RPK_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
ncol = length(RNA)
collist1 = c(1,seq(2,ncol-2, by=2),ncol-1,ncol)
collist2 = c(1,seq(3,ncol-2, by=2),ncol-1,ncol)
RNA1 <- RNA[,..collist1]
RNA2 <- RNA[,..collist2]
} else {

# For TL10/TL12
RNA1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/RPK/TL10_',MapMethod,'_featureCounts_multi_RPK_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNA2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/RPK/TL12_',MapMethod,'_featureCounts_multi_RPK_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
}


# Use RNA read count table without ncRNA for normalizing
if (Expt == 'TL14') {

RNAfornorm <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/RPK/',Expt,'_',MapMethod,'_featureCounts_multi_reads_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNAfornorm1 <- RNAfornorm[,..collist1]
RNAfornorm2 <- RNAfornorm[,..collist2]
} else {

# For TL10/TL12
RNAfornorm1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/RPK/TL10_',MapMethod,'_featureCounts_multi_reads_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNAfornorm2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/RPK/TL12_',MapMethod,'_featureCounts_multi_reads_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))
}


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
tokeep <- c(1,grep('_0m', samplenames)+1)
DT_tmp <- merge(HLs, RNA_RPKMnuc1[,..tokeep], by.x='Gene', by.y='GeneName') 
DT <- merge(DT_tmp, RNA_RPKMnuc2[,..tokeep], by.x='Gene', by.y='GeneName') 
# Get abundance fold change
# DT[, RNA_FC1 := log(TL10_0m_LRP/TL10_0m_WT, 2)]
# DT[, RNA_FC2 := log(TL12_0m_LRP/TL12_0m_WT, 2)]
colnames <- colnames(DT)
DT[, RNA_FC1 := log(get(colnames[11])/get(colnames[10]), 2)]
DT[, RNA_FC2 := log(get(colnames[13])/get(colnames[12]), 2)]



heavyStrand = c('MT-RNR2', 'MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB')
lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1', '7S')


# Heavy strand

DTsub <- DT[DT$Gene %in% heavyStrand]


# First just plot mut vs wt
colnames <- colnames(DTsub)
xsamps = colnames[grep('_0m_WT', colnames)]
ysamps = colnames[grep('_0m_LRP', colnames)]

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
plot='Abundance_LRPvsWT'
if (Expt == 'TL14') {
outpath = paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/Comparisons/')
extrainfo = paste0('_', MapMethod,'_')
} else {
outpath = paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/Comparisons/')
extrainfo = paste0('_',MapMethod,'_WTmethod',WTmethod,'mutmethod',mutmethod,'_')
}

pdf(paste0(outpath,plot,'_', MapMethod,'_heavy.pdf'), width = 3, height = 3.5)

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

# Save plot

plot = 'LRPeffect_vs_HL'

pdf(paste0(outpath,plot,extrainfo,'heavy_noHLcorrection.pdf'), width = 3, height = 3.5, useDingbats=FALSE) 


# Now mut/wt FC vs wt HL
xsamps = c('HL_wt1', 'HL_wt2')
ysamps = c('HL_FC1', 'HL_FC2')

xName = 'WT half-life'
yName = 'Half-life log2 LRPPRC KO / WT'

# For plotting WT/KO
ysamps = c('HL_recFC1', 'HL_recFC2')
yName = 'Half-life log2 LRPPRC WT / KO'

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

colptsize = 1


# LRP abundance vs WT
ylimits=c(-3, 1) # c(-4.1, -.25)
ylimits=c(-1,3)
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

text(min(xlimits)+(max(xlimits)-min(xlimits))/4, .5, labels = mylabel, cex = .8)

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
ysd <- apply(DTsub[,..ysamps], 1, sd)    

# Get ranges
xrange = rowMaxs(as.matrix(DTsub[,..xsamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..xsamps]), value=TRUE)

yrange = rowMaxs(as.matrix(DTsub[,..ysamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..ysamps]), value=TRUE)



# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
plot='AbundanceFC_vs_HL_FC'
pdf(paste0(outpath,plot,extrainfo,'heavy_noHLcorrection.pdf'), width = 3, height = 3.5)


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
plot='AbundanceFC'
pdf(paste0(outpath,plot,'_',MapMethod,'_barplot.pdf'), width = 3, height = 3.5)


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
ysd <- apply(DTsub[,..ysamps], 1, sd)    



# Get ranges
xrange = rowMaxs(as.matrix(DTsub[,..xsamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..xsamps]), value=TRUE)

yrange = rowMaxs(as.matrix(DTsub[,..ysamps]), value=TRUE) - rowMins(as.matrix(DTsub[,..ysamps]), value=TRUE)


# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
plot='AbundanceFC_vs_HL'
pdf(paste0(outpath,plot,extrainfo,'FC_light_noHLcorrection.pdf'), width = 3, height = 3.5)

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




# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/Correlations_foldChange_HL_RNA_LRPPRC.R')

