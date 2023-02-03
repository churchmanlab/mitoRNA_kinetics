library(data.table)
library('scales')


wid = 2 # width of plot 2.4
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2



NucCol = "Average" # "Fraction_in_Ribo_ave" 
MitoCols = c('TL8', 'TL9','TL11')

NucDT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/FractionRiboAssociatedNucK562_RIcheck.txt'),header=TRUE, stringsAsFactors=FALSE))
# FractionRiboAssociatedNucK562.txt used first, from Erik:  "I asked Robert to double check my "fractions new in the ribosome" for the cyto side (it's a bit more complicated as there is three states). He corrected the way I calculated it."

MitoDT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Fraction associated with the ribosomeMito.txt'), header=TRUE, stringsAsFactors=FALSE))

meanNucVals = NucDT[[NucCol]]
meanMitoVals = rowMeans(MitoDT[,..MitoCols])


# Check correlation between nuclear replicates
# pearson_r = cor.test(~r1NucVals+r2NucVals, method = c('pearson')) 
# r2 = (pearson_r$estimate)^2
# mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
# 
# pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/K562_OXPHOS_rep2vsrep1_', type,'.pdf'), width = 3, height = 3.5)
# 
# limits=c(0,max(x,y))
# plot(r1NucVals,r2NucVals, log='xy',cex.axis = 0.7,  pch = 16, xlab = paste0(NucCol, ' rep1'), ylab = paste0(NucCol, ' rep2')) #, ylim=limits, xlim=limits)
# # Diagonal line
# abline(0,1, lwd = 0.5, col='black')
# # Trendline
# #add linear trend
# abline(lm(r2NucVals~r1NucVals),col='black', lwd = 0.5, lty=2)
# 
# text(max(limits)/7, max(limits)*(4/5), labels = mylabel, cex = .8)
# 
# dev.off()

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/MitoNucFracRiboAssocBoxplot.pdf'), width = wid, height = ht)

xx=boxplot(meanMitoVals, meanNucVals, ylab = 'Fraction ribosome-associated',  xlab = '', names=c('Mito', 'NucOXPHOS'), pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt,outline=FALSE, las=las, ylim=c(0,1))


dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots_fracAssociated.r')

