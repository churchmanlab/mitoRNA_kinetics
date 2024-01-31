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

types=c("Half_life_whole_cell", "Half_life_poly_entry", "Half_life_nuc") 

meanNucVals = c()

for (type in types) {
NucCol = paste0(type,".MAP")


Nuc1DT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Brendan_K562_rates_rep1.csv'), sep=",",header=TRUE, stringsAsFactors=FALSE))

Nuc2DT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Brendan_K562_rates_rep2.csv'), sep=",",header=TRUE, stringsAsFactors=FALSE))



# Combine columns of interest for both replicates
coi=c(NucCol)
NucDT <- cbind(r1=Nuc1DT, r2=Nuc2DT[, ..coi])
newcoi=paste0(c('r1.','r2.'), NucCol)
NucDT <- NucDT[get(newcoi[1]) < 5000 & get(newcoi[2]) <5000 & get(newcoi[1]) != 'NA' & get(newcoi[2]) != 'NA']



nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

# Keep only genes of interest
NucgoiDT = NucDT[r1.Symbol %in% nucgenes]
# `%notin%` <- Negate(`%in%`)
# NucNonOXPHOS = NucDT[r1.Symbol %notin% nucgenes]

Meanr1NucVal=mean(NucgoiDT[[paste0('r1.',NucCol)]])
Meanr2NucVal=mean(NucgoiDT[[paste0('r2.',NucCol)]])
meanNucVal=mean(Meanr1NucVal, Meanr2NucVal)


meanNucVals=c(meanNucVals, meanNucVal)
}


# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/GetOXPHOSAverages.r')

