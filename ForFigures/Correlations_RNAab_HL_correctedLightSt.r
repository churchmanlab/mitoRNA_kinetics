library('scales')
library(data.table)
library(rlist)
library(Rfast)
library(matrixStats)

HLmeasure='HL_C_AverageAll_norm'
SDmeasure='HL_C_STDV_norm'
RNAabMeasure='AllRPK_norm' # TConlyRPK_norm TCunlabRPK_norm AllRPK_norm


HLs <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/Corrected_Light-Strand_Half-lives_Averages.csv', sep=',',header=TRUE, stringsAsFactors=FALSE))

RNA_TL1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL1_2020_02/RPK/TL1_MTnorRNA_t5MMinformed4_withDups_frag_all_TL1_240m_RPK.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

RNA_TL3 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/RPK/TL3_MTnorRNA_t5MTMMinformed6_withDups_frag_all_TL3_240m_RPK.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

RNA_TL4 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/RPK/TL4_MTnorRNA_t5MTMMinformed6_withDups_frag_all_TL4_240m_RPK.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

# Averages for RNA abundance
RNA_aves <- rowMeans(cbind(RNA_TL1[[RNAabMeasure]], RNA_TL3[[RNAabMeasure]], RNA_TL4[[RNAabMeasure]]))
# Std dev for RNA abundance
RNAsds <- rowSds(cbind(RNA_TL1[[RNAabMeasure]], RNA_TL3[[RNAabMeasure]], RNA_TL4[[RNAabMeasure]]))
# Make RNA data table
RNA <- data.table(Gene=paste0('MT-',RNA_TL1$GeneName), RNAave=RNA_aves, RNAsd=RNAsds)

# Merge RNAab with HL data
HLcois=c('Gene',HLmeasure, SDmeasure)
DT <- merge(HLs[, ..HLcois], RNA, by='Gene')

genes=DT$Gene


xName = 'Corrected half-life'
yName = 'antiND1-normalized abundance'

# Get x and y
x = DT[[2]]
y = DT[[4]]


# Get standard deviations
xsd <- DT[[3]]   
ysd <- DT[[5]]  



# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# Save plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/LightStrand_',RNAabMeasure,'_vs_', HLmeasure,'.pdf'), width = 3, height = 3.5)


pointsize = 1.2

genelist=DT$Gene
# 7S [1] "MT-ND6"            "MT-antiATP8-6-CO3" "MT-antiCO1"       
#  [4] "MT-antiCO2"        "MT-antiCYB"        "MT-antiND1"       
#  [7] "MT-antiND2"        "MT-antiND3"        "MT-antiND4L-4"    
# [10] "MT-antiND5"       

colors = c('white',rep('lightskyblue', 1), 'grey80',rep('pink1', 2), 'limegreen', rep('lightskyblue', 5))

xlimits=c(0,76)
ylimits=c(0,10)

# Half lives vs abundance
plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, xlim=xlimits, ylim=ylimits, main = RNAabMeasure)
# Diagonal line
# abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
abline(lm(y~x),col='black', lwd = 0.5, lty=2)

text(max(xlimits)/4, max(ylimits)*(4/5), labels = mylabel, cex = .8)

# Colors
for (i in c(1:length(genes))) {
points(x[which(genes==genelist[i])], y[which(genes==genelist[i])], pch = 21, cex = pointsize, bg = colors[i], col='black', lwd=.8)
}
# error bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_RNAab_HL_correctedLightSt.R')

