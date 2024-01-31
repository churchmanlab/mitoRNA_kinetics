library(data.table)
library('scales')


path = '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/'

filename = 'TLseqHalflives_corr_10000.txt'

# Read in data
DT <- data.table(read.table(paste0(path, filename), header=TRUE, sep = '|'))
# Remove first column which is NA
DT <- DT[, 2:ncol(DT)]
DT[, NS_July2021 := NULL]
# Rename first column
colnames(DT)[1] <- 'Txpt'
# Change first column type from integer to character (sapply(DT, typeof) to see each column type. And trim white space
DT$Txpt <- trimws(as.character(DT$Txpt))

# Reorder columns
# cols = colnames(DT)
# setcolorder(DT, cols[c(1,3, )]

# Order txpts in table
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
ord_index = append(ord_index, which(grepl(paste0('^',txpt,'$'), DT$Txpt)))
}
setorder(DT[, .r := order(ord_index)], .r)[, .r := NULL]

# xsamples <- c('ND1.normalized.half.life')
# ysamples <- c('ND1.normalized.Abundance')
#  
# xlimits = NULL # NULL c(100, 50000)
# ylimits = NULL # c(1000, 2000000)
# trans = 1
ptsize = 1.5
# colptsize = 1
# 
# txpt = DT$Txpt

pdf(paste0(path, 'HalfLifeCorrelationMatrix.pdf'), 
     width=5.5, # width=4.5, width=2, width=2.5
     height=5.5, #4
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
# pairs(DT[,2:ncol(DT)], lower.panel = DataPanel, upper.panel = NULL, gap=.4, labels=c('Low 4sU\nTL rep1', 'Low 4sU\nTL rep2','High 4sU\nTL rep1','High 4sU\nNS average', 'NS July2021','EtBr chase\nChujo et al'))
pairs(DT[,2:ncol(DT)], lower.panel = DataPanel, upper.panel = NULL, gap=.4, labels=c('Low 4sU\nTL rep1', 'Low 4sU\nTL rep2','High 4sU\nTL rep1','High 4sU\nNS average','EtBr chase\nChujo et al'))



dev.off()

# With MT-RNR2
# Get correlation
# pearson_r = cor.test(~x+y, method = c('pearson'))
# mylabel = bquote(R^2 == .(format(pearson_r$estimate^2, digits = 3)))
# # Plot
# plot(x,y, xlim = xlimits, ylim = ylimits, cex.axis = 0.7, main = 'All txpts', col = colors, pch = 16, cex = ptsize, xlab = xsamples[i], ylab = ysamples[i])
# # legend
# legend("bottomright", legend=c("RNR2", "Complex I", "Complex III", "Complex IV", "Complex V"), cex=.6, bty="n", text.col =  c(col1,col2,col5,col3,col4),  ncol=1)
# # Rsquared value
# text(x=diff(range(x))*.3+min(x),y=diff(range(y))*.9+min(y), labels = mylabel, cex = 1)
# # regression line
# abline(lm(y ~ x), lty = 2, lwd = .5)
# 
# # No log, no MT-RNR2
# x = x[2:length(x)]
# y = y[2:length(y)]
# txpt = txpt[2:length(txpt)]
# colors = colors[2:length(colors)]
# # Redo correlation without RNR2
# pearson_r = cor.test(~x+y, method = c('pearson'))
# mylabel = bquote(R^2 == .(format(pearson_r$estimate^2, digits = 3)))
# # Plot
# pairs[DT[,2:length(DT)],xlim = xlimits, ylim = ylimits, cex.axis = 0.7, main = 'No RNR2', col = colors, pch = 16, cex = ptsize, xlab = xsamples[i], ylab = ysamples[i])
# 
# #plot(x,y, xlim = xlimits, ylim = ylimits, cex.axis = 0.7, main = 'No RNR2', col = colors, pch = 16, cex = ptsize, xlab = xsamples[i], ylab = ysamples[i])
# # Legend
# legend("bottomright", legend=c("Complex I", "Complex III", "Complex IV", "Complex V"), cex=.6, bty="n", text.col =  c(col2,col5,col3,col4),  ncol=1)
# # Rsquared
# text(x=diff(range(x))*.3+min(x),y=diff(range(y))*.9+min(y), labels = mylabel, cex = 1)
# # regression line
# abline(lm(y ~ x), lty = 2, lwd = .5)

# 10 RPK line
# abline(h=10, lwd = 0.5, col='gray60', lty=2)
# abline(v=10, lwd = 0.5, col='gray60', lty=2)
## center line:
# if (type == 'RPKM' & set == 'all'){
# abline(0,1, lwd=0.6, col = 'gray40')
# 
# ## 2-fold lines:
# abline(log(2,10),1, lwd = 0.5, col='gray60', lty=2)
# abline(-log(2,10),1, lwd = 0.5, col='gray60', lty=2)
# }
# # abline(h=1000)
# # text(x = rlabelx, y = rlabely, labels = mylabel1, cex = .8)
# # Color points in complexes
# 
# # Set colors
# # colors = c('darkgoldenrod1', 'dodgerblue', 'indianred3', 'forestgreen')
# 
# heavyStrand = c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8-6','MT-CO3','MT-ND3','MT-ND4L-4','MT-ND5', 'MT-CYB')
# lightStrand = c('MT-antiCYB','MT-ND6', 'MT-antiND5','MT-antiND4L-4','MT-antiND3', 'MT-antiATP8-6-CO3', 'MT-antiCO2', 'MT-antiCO1', 'MT-antiND2', 'MT-antiND1')
# 
# if (set == 'all') {
# # Mito genes
# points(x[genes %in% heavyStrand], y[genes %in% heavyStrand],pch = 16, cex = colptsize, col = 'orange')
# points(x[genes %in% lightStrand], y[genes %in% lightStrand],pch = 17, cex = colptsize, col = 'red')
# # Mito rRNA
# points(x[grep('MT-RNR2', genes)], y[grep('MT-RNR2', genes)],pch = 16, cex = colptsize, col = 'purple')
# # Cyto rRNA
# points(x[grep('LSU_rRNA', genes)], y[grep('LSU_rRNA', genes)],pch = 16, cex = colptsize, col = 'green')
# points(x[grep('SSU_rRNA_18S', genes)], y[grep('SSU_rRNA_18S', genes)],pch = 16, cex = colptsize, col = 'dodgerblue')
# legend("bottomright", legend=c("heavy strand", "light strand", "Mito RNR2", "LSU rRNA", "SSU rRNA"), cex=.6, bty="n", text.col =  c('orange', 'red','purple', 'green', 'dodgerblue'),  ncol=1)
# 
# }
# 
# if (set == 'mito') {
# points(x, y,pch = c(rep(16, 10), rep(17, 10)), cex = 1, col = colors)
# legend("bottomright", legend=c('CI', 'CIII', 'CIV', 'CV'), cex=1, bty="n",text.col = c('darkgoldenrod1', 'dodgerblue', 'indianred', 'forestgreen'))
# text(x = .01*max(x), y = .8*max(y), labels = mylabel, cex = 1)
# }



# List genes increased in mito enrichment
# MitoEnriched = genes[y>1000 & y>x]
# write.table(MitoEnriched, file=paste0(path, 'Correlations/', yfile, '_vs_', xfile, '_', i,'.txt'),row.names=FALSE, sep=("\t"), quote=FALSE)


# identify(x, y, labels = genes)

# dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlotting/CorrelationPlot_HalfLifeComparison_matrix.R')

