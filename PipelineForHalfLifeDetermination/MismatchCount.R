

library(data.table)
library(rlist)
library('scales')
library('plotrix')

Folder = "HEK_TL13_LRPPRC_mitoriboIP_2022_08" 
set= 'MTnorRNA' # MT Nuc
Exp = 'TL13_B'
# Samples = c("0m_WT_tot_A","30m_WT_tot_A","60m_WT_tot_A","0m_LRP_tot_A","30m_LRP_tot_A","60m_LRP_tot_A","0m_WT_IP_A","30m_WT_IP_A","60m_WT_IP_A","0m_LRP_IP_A","30m_LRP_IP_A","60m_LRP_IP_A") # TL13_A
Samples = c("0m_WT_tot_B","30m_WT_tot_B","60m_WT_tot_B","0m_LRP_tot_B","30m_LRP_tot_B","60m_LRP_tot_B","0m_WT_IP_B","30m_WT_IP_B","60m_WT_IP_B","0m_LRP_IP_B","30m_LRP_IP_B","60m_LRP_IP_B") # TL13_B
# Samples = c("0m","90m_4sU","90m_6sG","180m_4sU6sG")
# c("0m","90m_4sU","90m_6sG_200","90m_6sG_1000","180m_4sU6sG_200","180m_4sU6sG_1000") # TL7
# c("0m_WT","0m_LRP","30m_WT","30m_LRP","60m_WT","60m_LRP","120m_WT","120m_LRP")
# c("0m_WT","0m_LRP","30m_WT","30m_LRP","60m_WT","60m_LRP","120m_WT","120m_LRP")
#c("0m_tot_A", "30m_tot_A", "60m_tot_A", "0m_IP_A", "30m_IP_A", "60m_IP_A", "0m_tot_B", "30m_tot_B", "60m_tot_B", "0m_IP_B", "30m_IP_B", "60m_IP_B") # c('DMSO', '1_NHC', '10_NHC', '100_NHC', '100_Molnu')
# c("0m_tot","30m_tot","60m_tot","0m_IP","30m_IP","60m_IP")
MapMethod = 't5MTMMinformed6' # t5MTMMinformed6 MMinformed4

path = paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',Folder,'/MMfrequency/')

DT <- data.table(read.table(paste0(path, Exp,'_', set, '_', MapMethod, '_withDups_mismatches.count'), sep = '\t', header=TRUE))

MMs = DT[[1]]

# cols = c('dodgerblue', 'gold', 'orange', 'red')
# cols = c('dodgerblue', 'forestgreen', 'gold', 'orange', 'red')
# cols = rep(c('darkblue','dodgerblue','forestgreen','lightgreen', 'gold', 'orange', 'red', 'brown'),2)
cols = rep(c('darkblue','dodgerblue','forestgreen', 'gold', 'orange', 'red'),2)
samplenames = colnames(DT)[2:length(DT)]

# Normalize to CG (least frequent mismatch across samples)
DT_norm <- copy(DT)
normcols <- list()
for (i in c(2:length(DT))) {
	normcols <- list.append(normcols, DT[[i]]/DT[[i]][2])
	}
DT_norm[, (samplenames) := normcols]



# Plot
pdf(paste0(path,Exp,'_', set,'_',MapMethod, 'AllMismatchesCount.pdf'), width = 12, height = 10)
par(mfrow=c(2,1), cex.lab = 1)

# Unnormalized
toPlot <- as.matrix(DT[,2:length(DT)])
barplot(t(toPlot), beside=TRUE,names.arg=MMs, cex.names=.8, xlab = 'Genome nt: Read nt', ylab = 'Number mismatches', main=paste0(set, ' ', MapMethod, ' No norm'), col=cols,legend.text=Samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8))

# Normalized
toPlot <- as.matrix(DT_norm[,2:length(DT)])
barplot(t(toPlot), beside=TRUE,names.arg=MMs, cex.names=.8, xlab = 'Genome nt: Read nt', ylab = 'Number mismatches norm to CG', main=paste0(set, ' ', MapMethod, ' Normalized to CG'), col=cols,legend.text=Samples, args.legend = list(x='topleft', bty = 'n', fill = cols, cex = .8))

dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/MismatchCount.R')

