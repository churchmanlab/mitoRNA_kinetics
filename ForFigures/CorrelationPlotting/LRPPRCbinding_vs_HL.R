

### Written 7/14/23 by M. Couvillion

### Get LRPPRC PAR-CLIP footprint info from https://www.nature.com/articles/s41467-017-01221-z#additional-information and correlate it with Erik's half-lives

library('scales')
library(data.table)
library(rlist)
library(Rfast)

input='footprints'
org='human'


directory1 <- 'HEK_TL10_LRPPRC_2022_07' # HEK_TL14_LRPPRC_2023_05
directory2 <- 'HEK_TL12_LRPPRC_2022_08'
MapMethod <- 't5MTMMinformed6'
Expt <- 'TL12'

# WT1 <- 'TL14_WT_A' # TL10_WT
# WT2 <- 'TL14_WT_B' # TL12_WT
# mut1 <- 'TL14_LRP_A'
# mut2 <- 'TL14_LRP_A'

WT1 <- 'TL10_WT' # TL10_WT
WT2 <- 'TL12_WT' # TL12_WT
mut1 <- 'TL10_LRP'
mut2 <- 'TL12_LRP'

WTmethod= '_PcMTnorRNA_' # _PcMTnorRNA_ (This for TL10 and TL12) # _ (TL14)
WTmod='halflives' #  halflives (This for TL10 and TL12) _30minEdit_halflives_ 30minEdit2_halflives # halflives (TL14)

mutmethod='_PcBrute1_' # PcMTnorRNA PcBrute1 PcBrute2 # PcBrute1 uses the RNR2 fraction new calculated from the average half-lives. PcBrute2 directly uses the average RNR2 fractions new # PcBrute1 (TL10, TL12) # _ (TL14)
mutmod='halflives' # halflives (for TL10, TL12, TL14) 30minEdit_halflives 30minEdit2_halflives

filename <- paste0(WT1,'_MT_',MapMethod,'_modeAll',WTmethod,'FracNew_',WTmod,'_corr_100000min_v2.txt')
# TL14_WT_B_MT_t5MTMMinformed6_modeAll_FracNew_halflives_corr_10000min_v2.txt
HLwt1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt1', 'HLnorm')))

filename <- paste0(mut1,'_MT_',MapMethod,'_modeAll',mutmethod,'FracNew_',mutmod,'_corr_100000min_v2.txt')
HLmut1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut1', 'HLnorm')))

filename <- paste0(WT2,'_MT_',MapMethod,'_modeAll',WTmethod,'FracNew_',WTmod,'_corr_100000min_v2.txt')
# TL14_WT_B_MT_t5MTMMinformed6_modeAll_FracNew_halflives_corr_10000min_v2.txt
HLwt2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_wt2', 'HLnorm')))

filename <- paste0(mut2,'_MT_',MapMethod,'_modeAll',mutmethod,'FracNew_',mutmod,'_corr_100000min_v2.txt')
HLmut2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory2,'/HalfLife/',filename), sep='\t',skip=1, stringsAsFactors=FALSE, col.names=c('Gene', 'Fit','HL_mut2', 'HLnorm')))



# Merge because genes are not in the same order above
HLs_tmp1 <- merge(HLwt1[, c('Gene', 'HL_wt1')], HLmut1[, c('Gene', 'HL_mut1')], by='Gene')
HLs_tmp2 <- merge(HLs_tmp1[, c('Gene', 'HL_wt1', 'HL_mut1')], HLwt2[, c('Gene', 'HL_wt2')], by='Gene')
HLs <- merge(HLs_tmp2[, c('Gene', 'HL_wt1', 'HL_mut1', 'HL_wt2')], HLmut2[, c('Gene', 'HL_mut2')], by='Gene')
# Get fold change
HLs[, HL_FC1 := log(HL_mut1/HL_wt1, 2)]
HLs[, HL_FC2 := log(HL_mut2/HL_wt2, 2)]
HLs[, HL_recFC1 := log(HL_wt1/HL_mut1, 2)]
HLs[, HL_recFC2 := log(HL_wt2/HL_mut2, 2)]

# Remove RNR1
HLs <- HLs[Gene != 'MT-RNR1']



##### Get LRPPRC binding info #####

if (input == 'footprints') {
LRP_FPbed <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/LRPPRC_PAR-CLIP_2017/41467_2017_1221_MOESM3_ESM.csv', sep=',', header=TRUE, stringsAsFactors=FALSE))
LRP_FPbed$Start <- as.integer(LRP_FPbed$Start) 
LRP_FPbed$End <- as.integer(LRP_FPbed$End) 

mitobed <- data.table(read.table('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/mchrM.bed.txt',  header=FALSE))

percovs=c()
genes=c()
for (i in c(1:nrow(mitobed))) {
gene=as.vector(mitobed[i]$V4)
gstart=mitobed[i]$V2
gend=mitobed[i]$V3
length=gend-gstart
strand=mitobed[i]$V6
fps <- LRP_FPbed[Start >= gstart & End <= gend & Strand == strand]

if (nrow(fps) > 0) {
if (strand == '+') {fps[, diff := End-Start]}
if (strand == '-') {fps[, diff := Start-End]}
fpcov = sum(fps$diff)
percov = fpcov/length*100
} else {
percov = 0
}
genes=c(genes,gene)
percovs=c(percovs, percov)
}

BIND = data.table(Gene=genes, Value=percovs)
}

if (input == 'PARCLIP') {

path='/Users/Mary/Desktop/Data/TimelapseSeq/LRPPRC_PAR-CLIP_2017/bedGraphs/'
normSamp <- 'Control_Hela_PARCLIP'
Lib = 'LRPPRC_Hela_PARCLIP'

suffixP = paste0('_chrMnorRNA_cpm5p_Pall')
suffixM = paste0('_chrMnorRNA_cpm5p_Pall')

sampDTP <- data.table(read.table(paste0(path,Lib,suffixP, '.bedGraph'), row.names=3))
normDTP <- data.table(read.table(paste0(path,normSamp,suffixP, '.bedGraph'), row.names=3))
sampDTM <- data.table(read.table(paste0(path,Lib,suffixM, '.bedGraph'), row.names=3))
normDTM <- data.table(read.table(paste0(path,normSamp,suffixM, '.bedGraph'), row.names=3))


DTP <- cbind(sampDTP, normDTP[[3]]) # column of interest is 3 because the original third col is now the row names (see above)
DTM <- cbind(sampDTM, normDTM[[3]]) # column of interest is 3 because the original third col is now the row names (see above)
names(DTP) <- c('chr', 'pos0', 'samp', 'norm')
names(DTM) <- c('chr', 'pos0', 'samp', 'norm')

# Subtract normalized control from normalized sample
DTP[, diff :=samp-norm]
DTM[, diff :=samp-norm]
DTP$diff[DTP$diff < 0] <- 0
DTM$diff[DTM$diff < 0] <- 0

mitobed <- data.table(read.table('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/hChrM_merged.bed',  header=FALSE))


rpkms=c()
genes=c()
for (i in c(1:nrow(mitobed))) {
strand=mitobed[i]$V6
if (strand == '+') {DT <- DTP}
if (strand == '-') {DT <- DTM}
gene=as.vector(mitobed[i]$V4)
gstart=mitobed[i]$V2
gend=mitobed[i]$V3
length=gend-gstart
rpkm <- sum(DT[pos0 >= gstart & pos0 <= gend]$diff)/length*1000
genes=c(genes,gene)
rpkms=c(rpkms, rpkm)
}
BIND = data.table(Gene=genes, Value=rpkms)
}


###################################

DT = merge(HLs, BIND, by='Gene')

heavyStrand = c('MT-RNR2', 'MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB')
lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1', '7S')





DTsub <- DT[DT$Gene %in% heavyStrand]

# Save plot

plot = paste0('LRPbinding_',org,'_vs_HL')

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/Comparisons/',plot,'_heavy.pdf'), width = 3, height = 3.5, useDingbats=FALSE) 


# Now mut/wt FC vs wt HL
xsamps = c('HL_wt1', 'HL_wt2')
ysamps = c('Value', 'Value')

xName = 'WT half-life'
if (input == 'footprints') {
yName = 'Percent of CDS bound by LRPPRC'}
if (input == 'PARCLIP') {
yName = 'RPKM of LRPPRC PAR-CLIP'}


# Get x and y
x = rowMeans(DTsub[,..xsamps])
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# remove RNR2
x = x[1:10]
y = y[1:10]
genes = genes[1:10]

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..ysamps], 1, sd)    
xsd = xsd[1:10]
ysd = ysd[1:10]

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

colptsize = 1


# LRP abundance vs WT
ylimits=NULL #c(-3, 1) # c(-4.1, -.25)
xlimits=c(10, 140)
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=ylimits, xlim=xlimits)
 # , ylim=limits, xlim=limits)

# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)


# Diagonal line
# abline(0,1, lwd = 0.5, col='black', lty=2)
# Trendline
# Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
points(x[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)], y[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)],pch = 21, cex = colptsize, bg='dodgerblue',col = 'black')
points(x[grep('MT-CYB', genes)], y[grep('MT-CYB', genes)],pch = 21, cex = colptsize, bg = 'aquamarine2',col = 'black')
points(x[grep('MT-CO1|MT-CO2|MT-CO3', genes)], y[grep('MT-CO1|MT-CO2|MT-CO3', genes)],pch = 21, cex = colptsize, bg = 'pink3',col = 'black')
points(x[grep('MT-ATP8-6', genes)], y[grep('MT-ATP8-6', genes)],pch = 21, cex = colptsize, bg = 'grey60',col = 'black')

# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# 7S RNA
# points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 16, cex = colptsize, col = 'brown')

text((min(x)+max(x))/4, max(y), labels = mylabel, cex = .8)

legend("topright", legend=c("CI", "CIII", "CIV", "CV"), cex=.8, bty="n", text.col =  c('dodgerblue', 'aquamarine2', 'pink3','grey60'),  ncol=1) 

dev.off()



plot = paste0('LRPeffect_vs_LRPbinding_',org)

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',directory1,'/Comparisons/',plot,'_heavy.pdf'), width = 3, height = 3.5, useDingbats=FALSE) 


# Now mut/wt FC vs wt HL
xsamps = c('Value', 'Value')
ysamps = c('HL_recFC1', 'HL_recFC2')

if (input == 'footprints') {
xName = 'Percent of CDS bound by LRPPRC'}
if (input == 'PARCLIP') {
xName = 'RPKM of LRPPRC PAR-CLIP'}
yName = 'Half-life log2 LRPPRC WT / KO'


# Get x and y
x = rowMeans(DTsub[,..xsamps])
y = rowMeans(DTsub[,..ysamps])
genes = DTsub$Gene

# remove RNR2
x = x[1:10]
y = y[1:10]
genes = genes[1:10]

# Get standard deviations
xsd <- apply(DTsub[,..xsamps], 1, sd)    
ysd <- apply(DTsub[,..ysamps], 1, sd)    
xsd = xsd[1:10]
ysd = ysd[1:10]

# Calculate correlation coefficients
pearson_r = cor.test(~x+y, method = c('pearson')) 
r2 = (pearson_r$estimate)^2
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

colptsize = 1


# LRP abundance vs WT
ylimits=NULL #c(-3, 1) # c(-4.1, -.25)
xlimits=NULL
plot(x,y, cex.axis = 1, col = 'white', xlab = xName, ylab = yName, ylim=ylimits, xlim=xlimits)
 # , ylim=limits, xlim=limits)

# std dev bars
arrows(x-xsd/2, y, x+(xsd/2), y, length=0.05, angle=90, code=3, lwd=.2)
arrows(x, y-ysd/2, x, y+(ysd/2), length=0.05, angle=90, code=3, lwd=.2)


# Diagonal line
# abline(0,1, lwd = 0.5, col='black', lty=2)
# Trendline
# Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
points(x[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)], y[grep('MT-ND1|MT-ND2|MT-ND3|MT-ND5|MT-ND4L-4', genes)],pch = 21, cex = colptsize, bg='dodgerblue',col = 'black')
points(x[grep('MT-CYB', genes)], y[grep('MT-CYB', genes)],pch = 21, cex = colptsize, bg = 'aquamarine2',col = 'black')
points(x[grep('MT-CO1|MT-CO2|MT-CO3', genes)], y[grep('MT-CO1|MT-CO2|MT-CO3', genes)],pch = 21, cex = colptsize, bg = 'pink3',col = 'black')
points(x[grep('MT-ATP8-6', genes)], y[grep('MT-ATP8-6', genes)],pch = 21, cex = colptsize, bg = 'grey60',col = 'black')

# Mito rRNA
points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# 7S RNA
# points(x[grep('7S', genes)], y[grep('7S', genes)],pch = 16, cex = colptsize, col = 'brown')

text((min(x)+max(x))/3, max(y), labels = mylabel, cex = .8)

legend("topright", legend=c("CI", "CIII", "CIV", "CV"), cex=.8, bty="n", text.col =  c('dodgerblue', 'aquamarine2', 'pink3','grey60'),  ncol=1) 

dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/LRPPRCbinding_vs_HL.R')

