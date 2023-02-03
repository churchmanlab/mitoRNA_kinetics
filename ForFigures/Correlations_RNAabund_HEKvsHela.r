library('scales')
library(data.table)
library(rlist)
library(Rfast)


RNR2='yes'

RNA_Hela1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/RPK/TL3_t5MTMMinformed6_featureCounts_multi_RPKM_noPseudo.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

RNA_Hela2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/RPK/TL4_t5MTMMinformed6_featureCounts_multi_RPKM_noPseudo.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))


RNA_HEK1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/RPK/TL10_t5MTMMinformed6_featureCounts_multi_RPKM_noPseudo.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

RNA_HEK2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/RPK/TL12_t5MTMMinformed6_featureCounts_multi_RPKM_noPseudo.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))


# Combine
RNA_Hela <- data.table(GeneName=RNA_Hela1$GeneName, Hela1=RNA_Hela1$TL3_0m,  Hela2=RNA_Hela2$TL4_0m)
RNA_HEK <- data.table(GeneName=RNA_HEK1$GeneName, HEK1=RNA_HEK1$TL10_0m_WT,  HEK2=RNA_HEK2$TL12_0m_WT)

# Keep only mito genes
# RNA_Hela <- RNA_Hela[grep('^MT-|^7S', GeneName)]
# RNA_HEK <- RNA_HEK[grep('^MT-|^7S', GeneName)]
# # Remove tRNA, RNR1
# RNA_Hela <- RNA_Hela[!grep('^MT-T|RNR1|anti', GeneName)]
# RNA_HEK <- RNA_HEK[!grep('^MT-T|RNR1|anti', GeneName)]

# Normalize to ND1
RNA_Hela[, Hela1norm := Hela1/RNA_Hela[GeneName=='MT-ND1']$Hela1]
RNA_Hela[, Hela2norm := Hela2/RNA_Hela[GeneName=='MT-ND1']$Hela2]
RNA_HEK[, HEK1norm := HEK1/RNA_HEK[GeneName=='MT-ND1']$HEK1]
RNA_HEK[, HEK2norm := HEK2/RNA_HEK[GeneName=='MT-ND1']$HEK2]


# heavyStrand = c('MT-RNR2', 'MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB', '7S')
# lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1')
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6', 'MT-RNR2')
if (RNR2 == 'no') {
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6')
}

colors = c(rep('lightskyblue', 5), 'lightgreen', rep('pink1', 3), 'grey80', 'gold1')
pointsize = 1.2


# subset

RNA_Hela <- RNA_Hela[RNA_Hela$GeneName %in% genelist]
RNA_HEK <- RNA_HEK[RNA_HEK$GeneName %in% genelist]

xsamps = c('Hela1norm', 'Hela2norm')
ysamps = c('HEK1norm', 'HEK2norm')

# Get x and y
x = rowMeans(RNA_Hela[,..xsamps]) 
y = rowMeans(RNA_HEK[,..ysamps])
genes = RNA_Hela$GeneName


# Get standard deviations
xsd <- apply(RNA_Hela[,..xsamps], 1, sd)    
ysd <- apply(RNA_HEK[,..ysamps], 1, sd)    

# Get ranges
xrange = rowMaxs(as.matrix(RNA_Hela[,..xsamps]), value=TRUE) - rowMins(as.matrix(RNA_Hela[,..xsamps]), value=TRUE)

yrange = rowMaxs(as.matrix(RNA_HEK[,..ysamps]), value=TRUE) - rowMins(as.matrix(RNA_HEK[,..ysamps]), value=TRUE)


# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))


pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/HEKvsHela_RNAab_RNR2',RNR2,'_SDbars.pdf'), width = 3, height = 3.5)

xName = 'HeLa S3 ND1-normalized abundance'
yName = 'HEK293T ND1-normalized abundance'

if (RNR2 == 'no') {
limits=c(0,4.5) # c(0,6) c(0,1200)
} else {
limits=c(0,80) # c(0,6) c(0,1200)
}


plot(x,y, cex.axis = 0.9, col = 'white',  pch = 21, xlab = xName, ylab = yName, ylim=limits, xlim=limits)
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







# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_RNAabund_HEKvsHela.R')

