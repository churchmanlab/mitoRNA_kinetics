library('scales')
library(data.table)
library(rlist)
library(Rfast)

cells=c('HelaModel', 'HelaRep2','HEKwt','HEKLRPPRC','K562', 'HEKLRPPRCwt') 
xcells=c('HelaModel', 'HelaRep2','HEKwt','HEKLRPPRC','K562', 'HEKLRPPRC') 
ycells=c('HelaModel', 'HelaRep2','HEKwt','HEKLRPPRC','K562', 'HEKwt') 
xCol='RPF' # HelaModel_RPF HelaRep2_RPF HEKwt_RPF HEKLRPPRC_RPF K562_RPF
yCols=c('HL', 'ProtSynth') 


DT <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Gene_expression_model2.txt', sep="\t",skip=2,header=TRUE, stringsAsFactors=FALSE))


genes =DT$Gene


# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/GeneExpressionModel_HLvsRPF_ProtSynthvsRPF.pdf'), width = 11, height = 4.5)

par(mfcol = c(2, 6))

for (i in c(1:length(xcells))) {
for (yCol in yCols) {

# Get x and y
x = DT[[paste0(xcells[i], '_', xCol)]]
y = DT[[paste0(ycells[i], '_', yCol)]]





# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))




pointsize = 1.2
# genelist = c('MT-ATP6', 'MT-ATP8', 'MT-CO1',  'MT-CO2',  'MT-CO3',  'MT-CYB',  'MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4',  'MT-ND4L', 'MT-ND5',  'MT-ND6')
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6', 'MT-RNR2')
# colors = c('limegreen', 'limegreen', 'indianred2', 'indianred2', 'indianred2', 'dodgerblue', rep('gold1',7))
colors = c(rep('lightskyblue', 5), 'lightgreen', rep('pink1', 3), 'grey80', 'gold1')


plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = paste0(xcells[i], ' ',xCol), ylab = paste0(ycells[i], ' ',yCol), main=cells[i])
# Diagonal line
# abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
abline(lm(y~x),col='black', lwd = 0.5, lty=2)

text(max(x)/3, max(y)*(4.5/5), labels = mylabel, cex = .8)

# Colors
for (j in c(1:length(genes))) {
points(x[which(genes==genelist[j])], y[which(genes==genelist[j])], pch = 21, cex = pointsize, bg = colors[j], col='black', lwd=.8)
}

}
}

dev.off()


# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/Correlations_GeneExpressionModel.r')

