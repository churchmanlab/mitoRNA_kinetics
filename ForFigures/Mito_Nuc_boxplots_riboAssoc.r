library(data.table)
library('scales')

# Ribosome association half-life in mito vs nuc
# From synthesis to ribo-associated and just the ribo-associated since export to cytosol.
# For nuclear side use Brendan's table (K562)
# Half_life_chr.MAP 13
# Half_life_nuc.MAP 11
# Half_life_nucexp_from_chr.MAP 6.90E-05
# Half_life_cyto.MAP  8.9
# Half_life_poly_entry.MAP 2.9
# Half_life_whole_cell.MAP 19
# Half_life_nucexp_from_nucdeg.MAP 26
# Half_life_nucdeg.MAP 89


wid = 2.4 # width of plot
ht = 3.5 # height of plot
txt = .9 # axis text
las = 2

types=c("ExportToRiboAssoc", "SynthToRiboAssoc") # ExportToRiboAssoc SynthToRiboAssoc

meanNucValss = list()

for (type in types) {
NucCol1 = "Half_life_poly_entry.MAP" 
NucCol2 = "Half_life_nuc.MAP" 

if (type=="SynthToRiboAssoc") {
NucCol = "Half_life_poly_entry_fromExport.MAP"} 
if (type=="ExportToRiboAssoc") {
NucCol = "Half_life_poly_entry.MAP" 
}

MitoCol = "Average_HL_Ribo"

Nuc1DT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Brendan_K562_rates_rep1.csv'), sep=",",header=TRUE, stringsAsFactors=FALSE))

Nuc2DT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Brendan_K562_rates_rep2.csv'), sep=",",header=TRUE, stringsAsFactors=FALSE))

# Add column that is the sum of export and ribo association
Nuc1DT[, Half_life_poly_entry_fromExport.MAP := rowSums(.SD), .SDcols = c(NucCol1, NucCol2)]
Nuc2DT[, Half_life_poly_entry_fromExport.MAP := rowSums(.SD), .SDcols = c(NucCol1, NucCol2)]

# Combine columns of interest for both replicates
coi=c(NucCol)
# NucDT <- cbind(r1=Nuc1DT, r2=Nuc2DT[, ..coi])
NucDT <- cbind(r1=Nuc1DT, r2=Nuc2DT) # Added later, just for writing output
newcoi=paste0(c('r1.','r2.'), NucCol)
NucDT <- NucDT[get(newcoi[1]) < 5000 & get(newcoi[2]) <5000 & get(newcoi[1]) != 'NA' & get(newcoi[2]) != 'NA']


MitoDT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Ribo and processing rates2.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))


# Get genes of interest
mitogenes = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6')

nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

# Keep only genes of interest
MitogoiDT = MitoDT[Gene %in% mitogenes]
NucgoiDT = NucDT[r1.Symbol %in% nucgenes]

`%notin%` <- Negate(`%in%`)
# NucNonOXPHOS = NucDT[r1.Symbol %notin% nucgenes]

MitoVals=MitogoiDT[[MitoCol]]
r1NucVals=NucgoiDT[[paste0('r1.',NucCol)]]
r2NucVals=NucgoiDT[[paste0('r2.',NucCol)]]
meanNucVals=rowMeans(cbind(r1NucVals, r2NucVals))
# r1NucNonOXPHOSVals=NucNonOXPHOS[[paste0('r1.',NucCol)]]
# r2NucNonOXPHOSVals=NucNonOXPHOS[[paste0('r2.',NucCol)]]
# meanNucNonOXPHOSVals=rowMeans(cbind(r1NucNonOXPHOSVals, r2NucNonOXPHOSVals))
genesNuc=NucgoiDT$r1.Symbol
tab=cbind(genesNuc, meanNucVals)
write.table(tab, file=paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/NucOXPHOS_filtered_',type,'_mean.txt'), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)


meanNucValss[[length(meanNucValss) + 1]] <- meanNucVals
}

write.table(NucgoiDT, file=paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/NucOXPHOS_filtered_allVals.txt'), sep='\t', quote=FALSE, row.names=FALSE)
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

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoVsNuc/MitoNucRiboAssocBoxplot.pdf'), width = wid, height = ht)

xx=boxplot(MitoVals, meanNucValss[[1]], meanNucValss[[2]], ylab = 'Half-life ribo association (min)',  xlab = '', names=c('Mito', 'NucOXPHOS\nfromExport', 'NucOXPHOS\nfromSynth'), pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt, log='y',outline=FALSE, las=las)


dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Mito_Nuc_boxplots_riboAssoc.r')

