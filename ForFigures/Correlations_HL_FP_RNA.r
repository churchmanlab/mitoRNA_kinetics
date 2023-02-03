library('scales')
library(data.table)
library(rlist)
library(Rfast)

RNR2='yes'
count='whole' # whole 100ntCorr
log='yes'

HL1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_1592min.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

HL2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_1592min.txt', sep='\t',header=TRUE, stringsAsFactors=FALSE))

HLs <- merge(HL1[, c('Gene', 'Normalized.Half.Life')], HL2[, c('Gene', 'Normalized.Half.Life')], by='Gene')

# FP <- data.table(read.table('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/HelaFracs/hMitoRP/RPK/HelaFracs_featureCounts_allSizes5p_ignore1stLast_multi_noDups_mito_tpm_noPseudo_ncRNA.txt', header=TRUE, stringsAsFactors=FALSE))
# if (FPtype == 'IP'){
# FP <- data.table(read.table('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/HelaIP/hMitoRP/RPK/HelaIP_featureCounts_allSizesAsite_CDSignore1stlast_multi_noDups_mito_tpm_noPseudo_ncRNA.txt', header=TRUE, stringsAsFactors=FALSE))
# }

# Separate out ATP8-6 and ND4L-4 for comparing to FPs
# HL <- rbind(HLs, HLs[5,], HLs[8,])
# HL[5,2] <- 'MT-ATP8'
# HL[8,2] <- 'MT-ND4L'
# HL[11,2] <- 'MT-ATP6'
# HL[12,2] <- 'MT-ND4'

# We use the SEM2 (polyA-tailing), SEM7 (ligation) and SEM9 (ligation) experiments, update 8/9/22: actually he uses SEM2 and SEM3, maybe SEM4, which all have higher coverage and all done using polyA-tailing (also used for Robert's modeling)
# RNA1 <- data.table(read.table('/Users/Mary/Dropbox/mito_nanopore_seq/results_files/SEM2_gene_counts_3prime_100nt_window_vs_whole_gene.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
# RNA2 <- data.table(read.table('/Users/Mary/Dropbox/mito_nanopore_seq/results_files/SEM3_gene_counts_3prime_100nt_window_vs_whole_gene.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
# # Add normalized column
# RNA2[, counts_3prime_ND1_norm := counts_100nt_3prime/RNA2[gene_name=='MT-ND1']$counts_100nt_3prime]
# RNA3 <- data.table(read.table('/Users/Mary/Dropbox/mito_nanopore_seq/results_files/SEM4_gene_counts_3prime_100nt_window_vs_whole_gene.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
# # Add normalized column
# RNA3[, counts_3prime_ND1_norm := counts_100nt_3prime/RNA3[gene_name=='MT-ND1']$counts_100nt_3prime]

# Update 8/11/22 Newly corrected tables (from SEM2 and SEM3)
RNA1 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/AbundanceData/total_RNA_polyA_tailing_rep1_gene_counts_corrected.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNA1[, count_100nt_3prime_corrected_ND1norm := count_100nt_3prime_corrected/RNA1[name_gene=='MT-ND1']$count_100nt_3prime_corrected]
RNA1[, count_whole_corrected_ND1norm := count_whole_corrected/RNA1[name_gene=='MT-ND1']$count_whole_corrected]

RNA2 <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/AbundanceData/total_RNA_polyA_tailing_rep2_gene_counts_corrected.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
RNA2[, count_100nt_3prime_corrected_ND1norm := count_100nt_3prime_corrected/RNA2[name_gene=='MT-ND1']$count_100nt_3prime_corrected]
RNA2[, count_whole_corrected_ND1norm := count_whole_corrected/RNA2[name_gene=='MT-ND1']$count_whole_corrected]

# Merge 
RNAs <- merge(RNA1[, c('name_gene', 'count_100nt_3prime_corrected_ND1norm', 'count_whole_corrected_ND1norm')], RNA2[, c('name_gene', 'count_100nt_3prime_corrected_ND1norm', 'count_whole_corrected_ND1norm')], by='name_gene')

# Separate out ATP8-6 and ND4L-4 *only if comparing to FPs
# RNA <- rbind(RNAs, RNAs[1,], RNAs[9,])
# RNA[1,1] <- 'MT-ATP8'
# RNA[9,1] <- 'MT-ND4L'
# RNA[16,1] <- 'MT-ATP6'
# RNA[17,1] <- 'MT-ND4'

# Illumina RNA data
# RNA <- data.table(read.table('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/HelaErik/RNAseq/RPK/HelaErik_featureCounts_CDSignore1stlast_multi_noDups_mito_RPK_noPseudo_ncRNA.txt', header=TRUE, stringsAsFactors=FALSE))

# IP <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/mitoRiboIP_abundances.txt', header=TRUE, stringsAsFactors=FALSE))

# FP_HL <- merge(FP, HL, by.x='GeneName', by.y='X.1')
# FP_HL_RNA <- merge(FP_HL, RNA, by.x='GeneName', by.y='gene_name')

RNA_HL <- merge(RNAs, HLs, by.x='name_gene', by.y='Gene')

# Remove RNR1 and ND6, tRNAs, anti genes
RNA_HL <- RNA_HL[name_gene != 'MT-RNR1']
RNA_HL <- RNA_HL[name_gene != 'MT-ND6']
RNA_HL <- RNA_HL[!grep('MT-T', name_gene)]
RNA_HL <- RNA_HL[!grep('anti', name_gene)]

limits = c(0, 130)
if (RNR2 == 'no') {
# remove RNR2
RNA_HL <- RNA_HL[name_gene != 'MT-RNR2']
limits=c(0,8)
}
if (log=='yes') {
limits=c(1,130)
}
genes = as.character(RNA_HL$name_gene)

# FP_IP <- merge(FP, IP, by='GeneName')
# IPgenes = as.character(FP_IP$GeneName)

if (count == '100nt') {
xSamp = c('Normalized.Half.Life.x','Normalized.Half.Life.y')
ySamp = c('count_100nt_3prime_corrected_ND1norm.x', 'count_100nt_3prime_corrected_ND1norm.y')
} else if (count == 'whole') {
xSamp = c('Normalized.Half.Life.x','Normalized.Half.Life.y')
ySamp = c('count_whole_corrected_ND1norm.x', 'count_whole_corrected_ND1norm.y')
}

# FP_HL_RNA[,..ySamp2] <- as.numeric(as.matrix(FP_HL_RNA[,..ySamp2]))

xName = 'ND1-normalized half-life'
yName = 'ND1-normalized abundance'
# yName2 = 'Relative RNA abundance (ONT SEM2,7,9)'
# yName3 = 'Relative mitoRiboIP abundance'

# Get x and y
x = rowMeans(RNA_HL[,..xSamp])
y = rowMeans(RNA_HL[,..ySamp])
# y2 = rowMeans(FP_HL_RNA[,..ySamp2])
# 
# x3 = rowMeans(FP_IP[,..xSamp])
# y3 = rowMeans(FP_IP[,..ySamp3])


# Get standard deviations
xsd <- apply(RNA_HL[,..xSamp], 1, sd)    
ysd <- apply(RNA_HL[,..ySamp], 1, sd)    

# Get ranges
xrange = rowMaxs(as.matrix(RNA_HL[,..xSamp]), value=TRUE) - rowMins(as.matrix(RNA_HL[,..xSamp]), value=TRUE)

yrange = rowMaxs(as.matrix(RNA_HL[,..ySamp]), value=TRUE) - rowMins(as.matrix(RNA_HL[,..ySamp]), value=TRUE)

# y2range = rowMaxs(as.matrix(FP_HL_RNA[,..ySamp2]), value=TRUE) - rowMins(as.matrix(FP_HL_RNA[,..ySamp2]), value=TRUE)
# 
# x3range = rowMaxs(as.matrix(FP_IP[,..xSamp]), value=TRUE) - rowMins(as.matrix(FP_IP[,..xSamp]), value=TRUE)

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

# pearson_r1 = cor.test(~log(x,10)+log(y1,10), method = c('pearson')) 
# mylabel1 = bquote(italic(r) == .(format(pearson_r1$estimate, digits = 2)))
# # pearson_r2 = cor.test(~x+y2, method = c('pearson')) 
# pearson_r2 = cor.test(~log(x,10)+log(y2,10), method = c('pearson')) 
# mylabel2 = bquote(italic(r) == .(format(pearson_r2$estimate, digits = 2)))
# # pearson_r3 = cor.test(~x3+y3, method = c('pearson')) 
# pearson_r3 = cor.test(~log(x3,10)+log(y3,10), method = c('pearson')) 
# mylabel3 = bquote(italic(r) == .(format(pearson_r3$estimate, digits = 2)))

# Save plot
if (log == 'no') {
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Abundance',count,'_vs_HL_dt1592min_RNR2',RNR2,'_SDbars.pdf'), width = 3, height = 3.5)}
if (log == 'yes') {
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Abundance',count,'_vs_HL_dt1592min_RNR2',RNR2,'_SDbars_log.pdf'), width = 3, height = 3.5)}


pointsize = 1.2
# genelist = c('MT-ATP6', 'MT-ATP8', 'MT-CO1',  'MT-CO2',  'MT-CO3',  'MT-CYB',  'MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4',  'MT-ND4L', 'MT-ND5',  'MT-ND6')
genelist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6', 'MT-RNR2')
# colors = c('limegreen', 'limegreen', 'indianred2', 'indianred2', 'indianred2', 'dodgerblue', rep('gold1',7))
colors = c(rep('lightskyblue', 5), 'lightgreen', rep('pink1', 3), 'grey80', 'gold1')

# Half lives vs abundance
if (log == 'no') {
plot(x,y, cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, ylim=limits, xlim=limits)}
if (log == 'yes') {
plot(x,y, log='xy',cex.axis = 0.7, col = 'white',  pch = 21, xlab = xName, ylab = yName, ylim=limits, xlim=limits)}

# Diagonal line
abline(0,1, lwd = 0.5, col='black')
# Trendline
#add linear trend
# abline(lm(y~x),col='black', lwd = 0.5, lty=2)

text(max(limits)/4, max(limits)*(4/5), labels = mylabel, cex = .8)

# Colors
for (i in c(1:length(genes))) {
points(x[which(genes==genelist[i])], y[which(genes==genelist[i])], pch = 21, cex = pointsize, bg = colors[i], col='black', lwd=.8)
}
# Range bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Correlations_HL_FP_RNA.R')

