library(data.table)
library('scales')


DT1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/ND1normalizedRiboAbundance_ForSpread.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

DTrpf <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Ribo and processing rates3.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

DT <- merge(DT1, DTrpf, by='Gene')
####################################################################################
######################### Cumulative fraction plots for complexes ###################

abundanceMed=median(DT$ND1_normalized_abundnace.y)
riboMed=median(DT$AVERAGE_ND1_normalized_abundnace_in_ribosome.y)
rpfMed=median(DT$Average_RPF_tp10k_ND1_normalized)


# Fold difference from median
DT[, log2AbMed := log(ND1_normalized_abundnace.y/abundanceMed, 2)]
DT[, log2RiboMed := log(AVERAGE_ND1_normalized_abundnace_in_ribosome.y/riboMed, 2)]
DTrpf[, log2rpfMed := log(Average_RPF_tp10k_ND1_normalized/rpfMed, 2)]

Genes=DT$Gene

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CumulativeFraction.pdf'), width=3, height=2.5, pointsize=11)

plot(1, xlim=c(-2.5,1.5), ylim=c(0,1), pch=16, type='n', xlab='Synth rate relative to median (log2)', ylab='Cumulative fraction of Complexes', col='black', lwd=2, yaxt='n')
# Make shaded box at 2-fold
polygon(c(-.5,-.5,.5,.5), c(0,1,1,0), col = alpha('black', .1), border = NA)
# Draw lines 
x1=DT$log2AbMed
x1genes=Genes[order(x1)]
x1=sort(x1)
x1=c(x1[1], x1) # To make the stair start at the 0 line
x2=DT$log2RiboMed
x2genes=Genes[order(x2)]
x2=sort(x2)
x2=c(x2[1], x2) # To make the stair start at the 0 line
x3=DTrpf$log2rpfMed
x3genes=Genes[order(x3)]
x3=sort(x3)
x3=c(x3[1], x3) # To make the stair start at the 0 line

y=c(0:(length(x1)-1))
y=y/max(y)

points(x3,y, pch=16, type='s', col='red', lwd=2)
points(x1,y, pch=16, type='s', col='black', lwd=2)
points(x2,y, pch=16, type='s', col='forestgreen', lwd=2)

legend('bottomright', legend=c('Abundance', 'RiboAbundance', 'RPFabundance'), text.col=c('black','forestgreen', 'red'), bty='n', cex=.5)


axis(2, at=c(0,.5, 1), labels=c(0, .5, 1))

dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CumulativeFraction.R')
