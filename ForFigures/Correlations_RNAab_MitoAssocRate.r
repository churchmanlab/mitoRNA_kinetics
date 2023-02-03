library('scales')
library(data.table)
library(rlist)
library(Rfast)

RNR2='no'
count='whole' # whole 100ntCorr

# Abundance
# Update 8/11/22 Newly corrected tables (from SEM2 and SEM3)
RNA1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/AbundanceData/total_RNA_polyA_tailing_rep1_gene_counts_corrected.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNA1[, count_100nt_3prime_corrected_ND1norm := count_100nt_3prime_corrected/RNA1[name_gene=='MT-ND1']$count_100nt_3prime_corrected]
RNA1[, count_whole_corrected_ND1norm := count_whole_corrected/RNA1[name_gene=='MT-ND1']$count_whole_corrected]

RNA2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/AbundanceData/total_RNA_polyA_tailing_rep2_gene_counts_corrected.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNA2[, count_100nt_3prime_corrected_ND1norm := count_100nt_3prime_corrected/RNA2[name_gene=='MT-ND1']$count_100nt_3prime_corrected]
RNA2[, count_whole_corrected_ND1norm := count_whole_corrected/RNA2[name_gene=='MT-ND1']$count_whole_corrected]

# Merge 
RNAs <- merge(RNA1[, c('name_gene', 'count_100nt_3prime_corrected_ND1norm', 'count_whole_corrected_ND1norm')], RNA2[, c('name_gene', 'count_100nt_3prime_corrected_ND1norm', 'count_whole_corrected_ND1norm')], by='name_gene')


# Mito Association Rank
# already averaged
Rank <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/MitoAssocRank.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))


# Merge both
DT <- merge(RNAs, Rank, by.x='name_gene', by.y='GENE')

# limits = c(0, 130)
# if (RNR2 == 'no') {
# # remove RNR2
# RNA_HL <- RNA_HL[name_gene != 'MT-RNR2']
# limits=c(0,8)
# }
# 
genes = as.character(DT$name_gene)


xSamp = c('AVERAGE.RANK')

if (count == '100nt') {
ySamp = c('count_100nt_3prime_corrected_ND1norm.x', 'count_100nt_3prime_corrected_ND1norm.y')
} else if (count == 'whole') {
ySamp = c('count_whole_corrected_ND1norm.x', 'count_whole_corrected_ND1norm.y')
}


xName = 'ND1-normalized rank'
yName = 'ND1-normalized abundance'

# Get x and y
x = DT[[xSamp]]
y = rowMeans(DT[,..ySamp])

# Get standard deviations
xsd <- DT$STDEV   
ysd <- apply(DT[,..ySamp], 1, sd)    


# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))


# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Abundance',count,'_vs_MitoRiboAssocRank_SDbars.pdf'), width = 3, height = 3.5)

limits=c(0,12)
pointsize = 1.2
# genelist = c('MT-ATP6', 'MT-ATP8', 'MT-CO1',  'MT-CO2',  'MT-CO3',  'MT-CYB',  'MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4',  'MT-ND4L', 'MT-ND5',  'MT-ND6')
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6', 'MT-RNR2')
# colors = c('limegreen', 'limegreen', 'indianred2', 'indianred2', 'indianred2', 'dodgerblue', rep('gold1',7))
colors = c(rep('lightskyblue', 5), 'lightgreen', rep('pink1', 3), 'grey80', 'gold1')

# Half lives vs abundance
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

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_RNAab_MitoAssocRate.r')

