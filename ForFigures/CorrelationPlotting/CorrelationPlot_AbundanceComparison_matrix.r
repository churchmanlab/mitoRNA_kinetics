library(data.table)
library('scales')
library('rlist')

path = '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/'

filename = 'All_abundances_2.txt'
NoDupsfilename = 'Mito1_MTMMinformed6_featureCounts_multi_noDups_RPKM_noPseudo_MT.txt'

# Read in data
DTdups <- data.table(read.table(paste0(path, filename), header=TRUE, sep = '\t', fill=TRUE))
DTnoDups <- data.table(read.table(paste0(path, NoDupsfilename), header=TRUE, sep = '\t'))
# DT <- merge(DTdups, DTnoDups[,c('GeneName','Mito1_noDups')], by='GeneName', All=FALSE)
DT <- DTdups
# Keep only heavy strand non-rRNA
genes = c('MT-ND1', 'MT-ND2', 'MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5','MT-CYB')
DT <- DT[GeneName %in% genes]
# Leave just one of each technique
# DT <- DT[,c('GeneName','Mito1', 'Mito1_noDups', 'Tot1','TL2_0m_noC','Total_Nanostrings_2', 'Chujo_et_al','Mito_minION_3', 'Total_minION_3')]
# DT <- DT[,c('GeneName','Mito1', 'Tot1','TL2_0m_noC','Total_Nanostrings_2', 'Chujo_et_al','Mito_minION_3', 'Total_minION_3')]
DT <- DT[,c('GeneName','Mito1', 'Tot1','TL3_0m','Total_Nanostrings_2', 'Chujo_et_al','Rep_1_count_whole_corrected')]
# Normalize to ND1
sampcols = colnames(DT[,2:ncol(DT)])
ND1normDT <- copy(DT)
ND1norm <- list()
for (i in c(1:length(sampcols))) {
	ND1norm <- list.append(ND1norm, ND1normDT[[sampcols[i]]]/ND1normDT[GeneName=='MT-ND1'][[sampcols[i]]])
	}
ND1normDT[, (sampcols) := ND1norm]


# Make sure order is correct
to_ord = c('MT-ND1', 'MT-ND2', 'MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5','MT-CYB')
# Set colors
col1='darkgoldenrod1'
col2='dodgerblue'
col3='deeppink3'
col4='grey60'
col5='aquamarine2'
colors = c(rep(col2,2), rep(col3,2), rep(col4,1), rep(col3,1), rep(col2,3), rep(col5,1))

ord_index = c()
for (txpt in to_ord) {
ord_index = append(ord_index, which(grepl(paste0('^',txpt,'$'), ND1normDT$GeneName)))
}
setorder(ND1normDT[, .r := order(ord_index)], .r)[, .r := NULL]

ptsize = 1.5

pdf(paste0(path, 'AbundanceCorrelationMatrix_withTL3.pdf'), 
     width=7, # width=4.5, width=2, width=2.5
     height=7, #4
     pointsize=11)

# Correlation panel
# panel.cor <- function(x, y) {
#     usr <- par("usr"); on.exit(par(usr))
#     par(usr = c(0, 1, 0, 1))
# #     r <- round(cor(x, y), digits=2)
#     pearson_r = cor.test(~x+y, method = c('pearson'))
# #     txt <- bquote(R^2 == .(format(r^2, digits = 3)))
#     ptxt <- bquote(R^2 == .(format(pearson_r$estimate^2, digits = 3)))
#     cex.cor <- 0.8/strwidth(ptxt)
# #     text(0.5, 0.5, txt, cex = cex.cor * r^2)
#     text(0.5, 0.5, ptxt, cex = cex.cor * pearson_r$estimate^2)
# }

# Customize data panel
DataPanel <- function(x, y) {
  points(x,y, pch = 21, bg = alpha(colors, .5), cex = ptsize)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  pearson_r = cor.test(~x+y, method = c('pearson'))
  ptxt <- bquote(R^2==.(format(pearson_r$estimate^2, digits = 2)))
  cex.cor <- 0.8/strwidth(ptxt)
  text(0.35, 0.9, ptxt, cex = 1)
  }

# Customize diagonal panel
# text.panel<-function(x,y) {
# labels=c('Low 4sU\nTL rep1', 'Low 4sU\nTL rep2','High 4sU\nTL rep1','High 4sU\nNDS average','EtBr chase\nChujo et al')
# cex = ptsize
# text(0.5, 0.5, labels, cex = cex)
# }

# Create plots
pairs(ND1normDT[,2:(length(sampcols)+1)], lower.panel = DataPanel, upper.panel = NULL, gap=.4, labels=sampcols)


dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/CorrelationPlot_AbundanceComparison_matrix.R')

