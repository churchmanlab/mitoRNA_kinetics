library('scales')
library(data.table)
library(rlist)
library(Rfast)


RNR2='yes'

HL_Hela1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_Hela1', 'HLnorm')))

HL_Hela2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_Hela2', 'HLnorm')))


HL_HEK1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/HalfLife/TL10_WT_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_HEK1', 'HLnorm')))

HL_HEK2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/HalfLife/TL12_WT_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_10000min_v2.txt', sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_HEK2', 'HLnorm')))


# Merge because genes are not in the same order above
HLs_tmp1 <- merge(HL_Hela1[, c('Gene', 'HL_Hela1')], HL_Hela2[, c('Gene', 'HL_Hela2')], by='Gene')
HLs_tmp2 <- merge(HLs_tmp1[, c('Gene', 'HL_Hela1', 'HL_Hela2')], HL_HEK1[, c('Gene', 'HL_HEK1')], by='Gene')
HLs <- merge(HLs_tmp2[, c('Gene', 'HL_Hela1', 'HL_Hela2', 'HL_HEK1')], HL_HEK2[, c('Gene', 'HL_HEK2')], by='Gene')

# Remove RNR1
HLs <- HLs[Gene != 'MT-RNR1']




# heavyStrand = c('MT-RNR2', 'MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB', '7S')
# lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1')
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6', 'MT-RNR2')
if (RNR2 == 'no') {
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6')
}

colors = c(rep('lightskyblue', 5), 'lightgreen', rep('pink1', 3), 'grey80', 'gold1')
pointsize = 1.2


# subset

DTsub <- HLs[HLs$Gene %in% genelist]

xsamps = c('HL_Hela1', 'HL_Hela2')
ysamps = c('HL_HEK1', 'HL_HEK2')

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


pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/HEKvsHela_RNR2',RNR2,'_SDbars.pdf'), width = 3, height = 3.5)

xName = 'HeLa S3 half-life (min)'
yName = 'HEK293T half-life (min)'

if (RNR2 == 'no') {
limits=c(0,150) # c(0,6) c(0,1200)
} else {
limits=c(0,1200) # c(0,6) c(0,1200)
}


plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, ylim=limits, xlim=limits)
# Diagonal line
abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
abline(lm(y~x),col='black', lwd = 0.5, lty=2)

text(max(limits)/4, max(limits)*(4/5), labels = mylabel, cex = .8)

# Colors
for (i in c(1:length(genes))) {
points(x[which(genes==genelist[i])], y[which(genes==genelist[i])], pch = 21, cex = pointsize, bg = colors[i], col='black', lwd=.8)
}
# Range bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)


if (RNR2 == 'no') {
cols=c('lightskyblue',  'lightgreen','pink1', 'grey80')
legend('bottomright', legend = c('Complex I', 'Complex III', 'Complex IV', 'Complex V'), text.col = 'black', bty='n', pch = 21, pt.bg=cols, cex=.7)
} else {
cols=c('lightskyblue',  'lightgreen','pink1', 'grey80', 'gold1')
legend('bottomright', legend = c('Complex I', 'Complex III', 'Complex IV', 'Complex V', 'RNR2'), text.col = 'black', bty='n', pch = 21, pt.bg=cols, cex=.7)
}
dev.off()







# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/Correlations_HL_HEKvsHela.R')

