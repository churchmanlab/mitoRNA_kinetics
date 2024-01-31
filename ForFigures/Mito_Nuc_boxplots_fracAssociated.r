library(data.table)
library('scales')


wid = 2 # width of plot 2.4
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2


toUse = 'tot-nuc' # cyto
# From Erik: OK, so the data seem to be correct and it all comes down to the nuc OXPHOS genes being turned over very slowly So there is a large uncertainty in the cyto half-life and total half-life where the cyto-half lives are longer than expected based on rthe cyto half-lives. The question is which do we trust more? I calculated the % ribosome-associated mRNA using either the cyto as is or calculated the cyto as (total-nuc).  So you now have two columns for each replicate I tend to believe the longer total turnover half-lives more than the shorter cyto half-lives.

if (toUse == 'cyto') {
NucCol = "Average.Fraction in the polysome_Cyto" # "Fraction_in_Ribo_ave" 
}
if (toUse == 'tot-nuc') {
NucCol = "Average.Fraction in the polysome_Tot" # "Fraction_in_Ribo_ave" 
}

MitoCols = c('TL8', 'TL9','TL11')

NucDT <- fread(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Nuc_OXPHOS_Fraction in the polysome_update2024_0130.csv'),header=TRUE, stringsAsFactors=FALSE)
# FractionRiboAssociatedNucK562.txt used first, from Erik:  "I asked Robert to double check my "fractions new in the ribosome" for the cyto side (it's a bit more complicated as there is three states). He corrected the way I calculated it."
# Old data, bad gtf: /Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/FractionRiboAssociatedNucK562_RIcheck.txt
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


# Write table for Table S5
# the fraction ribosome-associated, the whole cell half-lives, nuclear half-lives, and polysome association rates
NucDT <- NucDT[order(Symbol)]
NucDT = NucDT[-c(1:2),] # First two rows are mean and median of column
write.table(NucDT[, c('Gene','Symbol','T.Fraction in the polysome_Tot', 'U.Fraction in the polysome_Tot', 'Average.Fraction in the polysome_Tot')], file=paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/NucOXPHOS_fractionInPolysome.txt'), sep='\t', quote=FALSE, row.names=FALSE, col.names=c('Gene','Gene symbol', 'Fraction in polysome Rep T', 'Fraction in polysome Rep U', 'Fraction in polysome average'))


file=paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/NucOXPHOS_filtered_',type,'_mean.txt'), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)



pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/MitoNucFracRiboAssocBoxplot_use',toUse,'.pdf'), width = wid, height = ht)

xx=boxplot(meanMitoVals, meanNucVals, ylab = 'Fraction ribosome-associated',  xlab = '', names=c('Mito', 'NucOXPHOS'), pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt,outline=FALSE, las=las, ylim=c(0,1))


dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots_fracAssociated.r')

