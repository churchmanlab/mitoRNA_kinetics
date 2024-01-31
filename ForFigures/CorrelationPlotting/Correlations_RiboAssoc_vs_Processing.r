library('scales')
library(data.table)
library(rlist)
library(Rfast)
library(matrixStats)

# Averages
RiboAssocRateName='Average_HL_Ribo'
ProcessingRateName='Procesing_Rates_min'

# Individual measurements
RiboAssocRateName1='TL8_HL_Ribo'
RiboAssocRateName2='TL9_HL_Ribo'
RiboAssocRateName3='TL11_HL_Ribo'


DT <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Ribo and processing rates.txt',header=TRUE, stringsAsFactors=FALSE))

genes=DT$Gene

# Averages for ribosome association rate
RiboAssocRate <- DT[[RiboAssocRateName]]
# Processing rate
ProcessingRate <- DT[[ProcessingRateName]]


# Specify axis labels
yName = 'Processing rate (min)'
xName = 'Mitoribosome association rate (min)'

# Get x and y
y = ProcessingRate
x = RiboAssocRate


# Get standard deviation for ribosome association rate
xsd <- rowSds(cbind(DT[[RiboAssocRateName1]], DT[[RiboAssocRateName2]], DT[[RiboAssocRateName3]]))


# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoriboAssocRate_vs_ProcessingRate.pdf'), width = 3, height = 3.5)


pointsize = 1.2

colors = c('grey80',rep('pink1', 3), 'limegreen', rep('lightskyblue',5))

# [1] "MT-ATP8-6" "MT-CO1"    "MT-CO2"    "MT-CO3"    "MT-CYB"    "MT-ND1"   
# [7] "MT-ND2"    "MT-ND3"    "MT-ND4L-4" "MT-ND5"   



ylimits=c(0,17)
xlimits=c(0,17)

# Half lives vs abundance
plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, xlim=xlimits, ylim=ylimits)
# Diagonal line
abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
# abline(lm(y~x),col='black', lwd = 0.5, lty=2)

# text(max(xlimits)/4, max(ylimits)*(4/5), labels = mylabel, cex = .8)

# Colors
points(x, y, pch = 21, cex = pointsize, bg = colors, col='black', lwd=.8)

# error bars
arrows(x-xsd/2, y, x+xsd/2, y, length=0.05, angle=90, code=3, lwd=.2)

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/Correlations_RiboAssoc_vs_Processing.R')

