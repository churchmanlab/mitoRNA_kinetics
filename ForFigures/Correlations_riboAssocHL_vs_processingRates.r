library('scales')
library(data.table)
library(rlist)
library(Rfast)
library(matrixStats)


end='5' # 5 3

# Get columns
xCol='Average_HL_Ribo'
xsdCol='STDV_HL_Ribo'
if (end == '5'){
yCol='Processing_HL_Average_5'
ysdCol='Processing_HL_STDV_5'
}
if (end == '3'){
yCol='Processing_HL_Average_3'
ysdCol='Processing_HL_STDV_3'
}



DT <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Ribo and processing rates3.txt',header=TRUE, stringsAsFactors=FALSE))

genes=DT$Gene

# Averages
y <- DT[[yCol]]
x <- DT[[xCol]]


# Specify axis labels
xName = 'Mitoribosome association half-life (min)'
yName = paste0(end,"' processing half-life (min)")

# Get sds
xsd=DT[[xsdCol]]
ysd=DT[[ysdCol]]

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/',end,'pProcessingHL_vs_riboAssocHL.pdf'), width = 3, height = 3.5)


pointsize = 1.2

colors = c('grey80',rep('pink1', 3), 'limegreen', rep('lightskyblue',5))

# [1] "MT-ATP8-6" "MT-CO1"    "MT-CO2"    "MT-CO3"    "MT-CYB"    "MT-ND1"   
# [7] "MT-ND2"    "MT-ND3"    "MT-ND4L-4" "MT-ND5"   


ylimits=c(0,25)
xlimits=c(0,25)

# Plot
plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, xlim=xlimits, ylim=ylimits, cex.lab=.8)
# Diagonal line
# abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
abline(lm(y~x),col='black', lwd = 0.5, lty=2)

text(max(xlimits)/4, max(ylimits)*(9.5/10), labels = mylabel, cex = .8)

# Colors
points(x, y, pch = 21, cex = pointsize, bg = colors, col='black', lwd=.8)

# error bars
arrows(x-xsd/2, y, x+xsd/2, y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+ysd/2, length=0.05, angle=90, code=3, lwd=.2)

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_riboAssocHL_vs_processingRates.R')

