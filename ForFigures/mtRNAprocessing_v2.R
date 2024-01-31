#!/n/app/R/4.0.1/bin/Rscript

# Plot processing data at 5' and 3' ends of mito transcripts

# Written by: M. Couvillion
# Date: March 29, 2022

# Use: sbatch -p short -t 0-00:30 -e logs/plotProcessing.err -o logs/plotProcessing_v2.log --mem=30G --wrap="/n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/mtRNAprocessing_v2.R LRPPRC_WT WT1_2 WT2_2 RPF"

# Use: sbatch -p short -t 0-00:30 -e logs/plotProcessing.err -o logs/plotProcessing.log --mem=30G --wrap="/n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/mtRNAprocessing_v2.R H_6 H_6_1 H_6_2 _chrMall_noDups RPF"
Exp <- 'H_6'
S1 <- 'H_5_1'
S1 <- 'H_6_1'
S2 <- 'H_6_2'
suffix <- '_chrMall_noDups'
type <- 'RPF'

# sbatch -p short -t 0-00:30 -e logs/plotProcessing.err -o logs/plotProcessing.log --mem=30G --wrap="/n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/mtRNAprocessing_v2.R Hela_IP Hela_IP_8U_1 Hela_IP_8U_2 _chrMall_noDups RPF"

# Use: sbatch -p short -t 0-00:30 -e logs/plotProcessing.err -o logs/plotProcessing.log --mem=30G --wrap="/n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/mtRNAprocessing_v2.R Fibro_RPF Fibro_1 Fibro_2 _chrMall_noDups RPF"

# sbatch -p short -t 0-00:30 -e logs/plotProcessing.err -o logs/plotProcessing.log --mem=30G --wrap="/n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/mtRNAprocessing_v2.R Fibro_RNA Fibro_1 Fibro_2 _Aligned.Mito.noDups RNA"

library(data.table)
library('scales')
library(Rfast)

args <- commandArgs(trailingOnly = TRUE)
# transcript = args[1]

# Read in parameters (data)

#path = '/n/groups/churchman/mc348/TimelapseSeq/TL8_mitoRiboIP_2021_12/' 
Exp <- args[1]
S1 <- args[2] # Hela_IP_8U Fibro H_6
S2 <- args[3]
suffix <- args[4] # _Aligned.Mito.noDups
type <- args[5] # RNA RPF

# Define ends for plus strand
txpts = c('ND1', 'ND2', 'CO1', 'CO2', 'ATP8_6', 'CO3', 'ND3', 'ND4L_4', 'ND5', 'CYB')
# strands = c(rep('+', 9), '-', '+')
fivepends = c(3305, 4470, 5901, 7586, 8365, 9207, 10059, 10470, 12337, 14747)
if (type == 'RPF') {
threepends = c(4262, 5511,7514, 8294, 9206, 9990, 10404, 12137, 14740,  15887) } # ND5 is problematic because there are no reads at the very 5' end of the transcript
if (type == 'RNA') {
threepends = c(4262, 5511,7514, 8294, 9206, 9990, 10404, 12137, 14745, 15887) 
}

Cols = c('name', 'samflag', 'chr', 'start', 'mapQ', 'cigar', 'seq', 'fiveP', 'threeP', 'pA', 'CCA', 'endinA', 'numMM')

sammP1 <- data.table(read.table(paste0(S1,suffix, '_P.samm'), col.names=Cols))
sammM1 <- data.table(read.table(paste0(S1, suffix, '_M.samm'), col.names=Cols))
sammP2 <- data.table(read.table(paste0(S2, suffix, '_P.samm'), col.names=Cols))
sammM2 <- data.table(read.table(paste0(S2, suffix, '_M.samm'), col.names=Cols))

# Make columns numeric
sammP1[, fiveP := as.numeric(fiveP)]
sammP1[, threeP := as.numeric(threeP)]
sammM1[, fiveP := as.numeric(fiveP)]
sammM1[, threeP := as.numeric(threeP)]

sammP2[, fiveP := as.numeric(fiveP)]
sammP2[, threeP := as.numeric(threeP)]
sammM2[, fiveP := as.numeric(fiveP)]
sammM2[, threeP := as.numeric(threeP)]


# Set leniency for being considered processed
a=2
# Set distances required to consider spanning junction
u=3
d=3

fivePunproc1=c()
threePunproc1=c()
fivePunproc2=c()
threePunproc2=c()
totals=c()
totals5p=c()
totals3p=c()
# And for keeping read counts to write table
fivePspan1s=c()
fivePspan2s=c()
fivePproc1s=c()
fivePproc2s=c()
threePspan1s=c()
threePspan2s=c()
threePproc1s=c()
threePproc2s=c()

for (i in c(1:length(txpts))) {
txpt = txpts[i]
# Find reads that span 5' junction
fivePspanDT1 = sammP1[fiveP < (fivepends[i]-u) & threeP > (fivepends[i]+d) & endinA == 'no']
# Find reads that are processed at 5' end
fivePprocDT1 = sammP1[fiveP >= (fivepends[i]-a) & fiveP <= (fivepends[i]+a)]
# Find reads that span 3' junction
threePspanDT1 = sammP1[fiveP < (threepends[i]-u) & threeP > (threepends[i]+d) & endinA == 'no']
# Find reads that are processed at 3' end
threePprocDT1 = sammP1[threeP >= (threepends[i]-a) & threeP <= (threepends[i]+a)]

# Rep 2
# Find reads that span 5' junction
fivePspanDT2 = sammP2[fiveP < (fivepends[i]-u) & threeP > (fivepends[i]+d) & endinA == 'no']
# Find redtads that are processed at 5' end
fivePprocDT2 = sammP2[fiveP >= (fivepends[i]-a) & fiveP <= (fivepends[i]+a)]
# Find reads that span 3' junction
threePspanDT2 = sammP2[fiveP < (threepends[i]-u) & threeP > (threepends[i]+d) & endinA == 'no']
# Find reads that are processed at 3' end
threePprocDT2 = sammP2[threeP >= (threepends[i]-a) & threeP <= (threepends[i]+a)]


# Get read counts
fivePspan1 = nrow(fivePspanDT1)
fivePproc1 = nrow(fivePprocDT1)
threePspan1 = nrow(threePspanDT1)
threePproc1 = nrow(threePprocDT1)
fivePspan2 = nrow(fivePspanDT2)
fivePproc2 = nrow(fivePprocDT2)
threePspan2 = nrow(threePspanDT2)
threePproc2 = nrow(threePprocDT2)
# Keep read counts
fivePspan1s=c(fivePspan1s, fivePspan1)
fivePspan2s=c(fivePspan2s, fivePspan2)
fivePproc1s=c(fivePproc1s, fivePproc1)
fivePproc2s=c(fivePproc2s, fivePproc2)
threePspan1s=c(threePspan1s, threePspan1)
threePspan2s=c(threePspan2s, threePspan2)
threePproc1s=c(threePproc1s, threePproc1)
threePproc2s=c(threePproc2s, threePproc2)

# Get totals
fiveP1=fivePspan1+fivePproc1
threeP1=threePspan1+threePproc1
fiveP2=fivePspan2+fivePproc2
threeP2=threePspan2+threePproc2
totals=c(totals, fiveP1, fiveP2, threeP1, threeP2)
totals5p=c(totals5p, fiveP1, fiveP2)
totals3p=c(totals3p, threeP1, threeP2)

# Get fraction unprocessed
fivePunprocFrac1 = fivePspan1/(fivePspan1 + fivePproc1)
threePunprocFrac1 = threePspan1/(threePspan1 + threePproc1)
fivePunprocFrac2 = fivePspan2/(fivePspan2 + fivePproc2)
threePunprocFrac2 = threePspan2/(threePspan2 + threePproc2)

# Store data
fivePunproc1 = c(fivePunproc1, fivePunprocFrac1)
threePunproc1 = c(threePunproc1, threePunprocFrac1)
fivePunproc2 = c(fivePunproc2, fivePunprocFrac2)
threePunproc2 = c(threePunproc2, threePunprocFrac2)
}



# Make barplot (5' and 3' together)
pdf(paste0(Exp,'_', type, '_fractionUnprocessed.pdf'), width = 8, height = 5.5)

ratios = matrix(c(fivePunproc1,fivePunproc2, threePunproc1,threePunproc2), nrow = 4, byrow=TRUE)
# counts = matrix(c(IP_pAcounts, IP_spancounts, tot_pAcounts, tot_spancounts), nrow = 4, byrow=TRUE)
# endcounts = matrix(c(IP_endcounts, tot_endcounts), nrow = 4, byrow=TRUE)


cols=c('black','black','red','red')

xx = barplot(ratios, beside=TRUE, names.arg=txpts, main=paste0(type, ' unprocessed ratios'), col=cols, cex.names = .8, ylab = 'Fraction unprocessed', ylim=c(0,1), legend.text=c("5'", "3'"), args.legend = list(x='topright', bty = 'n', fill = c('black','red'), cex = 1.2)) # 

# Add count as text to top of bars
text(x = as.vector(xx), y=0.8, label=totals, srt=90)
# text(x = xx, y = ratios+.07, label = counts, pos = 3, cex = 0.7, col = "black")

dev.off()

# For writing data to file
row.names(ratios) <- c('5p1', '5p2', '3p1', '3p2')
colnames(ratios) <- txpts
write.table(ratios, file=paste0(Exp,'_', type, '_fractionUnprocessed.txt'), quote=FALSE)




# Make barplot (5' and 3' separate)
pdf(paste0(Exp,'_', type, '_fractionUnprocessed_5p3pSeparate.pdf'), width = 7, height = 5.5) #  width = 8 for separate bars

ratios1 = matrix(c(fivePunproc1,fivePunproc2), nrow = 2, byrow=TRUE)
ratios2 = matrix(c(threePunproc1,threePunproc2), nrow = 2, byrow=TRUE)
ratiosave5p = (colMeans(ratios1))
ratiosave3p = (colMeans(ratios2))
maxs5p = (colMaxs(ratios1, value=TRUE))
mins5p = (colMins(ratios1, value=TRUE))
maxs3p = (colMaxs(ratios2, value=TRUE))
mins3p = (colMins(ratios2, value=TRUE))
ranges5p = maxs5p - mins5p
ranges3p = maxs3p - mins3p
# counts = matrix(c(IP_pAcounts, IP_spancounts, tot_pAcounts, tot_spancounts), nrow = 4, byrow=TRUE)
# endcounts = matrix(c(IP_endcounts, tot_endcounts), nrow = 4, byrow=TRUE)

cols=c('orange','orange')

## Reps plotted separately
# xx=barplot(ratios1, beside=TRUE, names.arg=txpts, main=paste0(type, " unprocessed ratios 5'"), col=cols, cex.names = .8, ylab = 'Fraction unprocessed', ylim=c(0,1)) # 
# # Add count as text to top of bars
# text(x = as.vector(xx), y=0.8, label=totals5p, srt=90)
# 
# xx=barplot(ratios2, beside=TRUE, names.arg=txpts, main=paste0(type, " unprocessed ratios 3'"), col=cols, cex.names = .8, ylab = 'Fraction unprocessed', ylim=c(0,1)) # 
# # Add count as text to top of bars
# text(x = as.vector(xx), y=0.8, label=totals3p, srt=90)

## Reps plotted together with range bars

xx=barplot(ratiosave5p, beside=TRUE, names.arg=txpts, main=paste0(type, " unprocessed ratios 5'"), col=cols, cex.names = .8, ylab = 'Fraction unprocessed', ylim=c(0,1)) # 
# range bars
arrows(xx, ratiosave5p-(ranges5p/2), xx, ratiosave5p+(ranges5p/2), length=0.03, angle=90, code=3, lwd=.2)

xx=barplot(ratiosave3p, beside=TRUE, names.arg=txpts, main=paste0(type, " unprocessed ratios 3'"), col=cols, cex.names = .8, ylab = 'Fraction unprocessed', ylim=c(0,1)) # 
arrows(xx, ratiosave3p-(ranges3p/2), xx, ratiosave3p+(ranges3p/2), length=0.03, angle=90, code=3, lwd=.2)

# Write counts for processed and unprocessed into table
ctTab = cbind(txpts, fivePspan1s, fivePproc1s, fivePspan2s, fivePproc2s, threePspan1s, threePproc1s, threePspan2s, threePproc2s)

colnames=c('Gene', paste0(Exp, c('_Span_5p_1',  '_Proc_5p_1', '_Span_5p_2', '_Proc_5p_2', '_Span_3p_1', '_Proc_3p_1', '_Span_3p_2', '_Proc_3p_2')))

write.table(ctTab, file=paste0(Exp,'_', type, '_processing_Counts.txt'), quote=FALSE, col.names=colnames)


dev.off()


