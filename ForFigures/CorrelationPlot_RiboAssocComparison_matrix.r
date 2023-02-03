library(data.table)
library('scales')


path = '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/'

fn = 'Combined_Fit_all_rates_20221122_twostate_Least-squares.txt'

# Read in data
DT <- data.table(read.table(paste0(path, fn), header=TRUE))


#          Gene       TL_9      TL_11         TL8
#  1:    MT-ND3 0.21703372 0.07546897 0.152553474
#  2:    MT-ND1 0.13630634 0.07024419 0.055113747
#  3:    MT-CYB 0.21871398 0.10248526 0.118774809
#  4:    MT-CO3 0.32568063 0.14566436 0.195859855
#  5:    MT-CO1 0.09901146 0.03659377 0.052356792
#  6: MT-ND4L-4 0.23295322 0.09178215 0.135657323
#  7:   MT-RNR2 0.01405373 0.00011884 0.003124897
#  8:    MT-ND5 0.26166126 0.13227003 0.127262487
#  9:    MT-CO2 0.24140788 0.09637931 0.124822102
# 10: MT-ATP8-6 0.23124263 0.09039087 0.163507841
# 11:    MT-ND2 0.16519850 0.09084312 0.093060813


# Reset column order
setcolorder(DT, c('Gene', 'TL8', 'TL_9', 'TL_11'))
sampnames = c('TL8', 'TL9', 'TL11')

# Set colors
col1='darkgoldenrod1' # 
col2='dodgerblue' # CI
col3='deeppink3' # CIV
col4='grey60' # CV
col5='aquamarine2' # CIII
colors = c(rep(col2,2), rep(col5,1), rep(col3,2), col2, col1, col2, col3, col4, col2)


ptsize = 1.5

pdf(paste0(path, 'RiboAssocCorrelationMatrix.pdf'), 
     width=3.5, # width=4.5, width=2, width=2.5
     height=3.5, #4
     pointsize=11)


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

# Create plots

pairs(DT[,2:ncol(DT)], lower.panel = DataPanel, upper.panel = NULL, gap=.4, labels=sampnames)


dev.off()



# Make barplot of the halflives (ln(2)/rate)

# Half lives
DT[,TL8_HL := log(2)/TL8]
DT[,TL9_HL := log(2)/TL_9]
DT[,TL11_HL := log(2)/TL_11]

ysamps = c('TL8_HL','TL9_HL', 'TL11_HL')

# Get x and y
y = rowMeans(DT[,..ysamps])
ysd <- apply(DT[,..ysamps], 1, sd)    

ord=order(y)
y=y[ord]
genesord=DT$Gene[ord]
ysdord = ysd[ord]

# Save plot
pdf(paste0(path, 'RiboAssoc_barplot.pdf'), width = 4, height = 3.5)

xx=barplot(y, log='y',ylim=c(1,260), names=genesord, las=2, cex.names=.7, ylab='Ribo unbound half-life (min)', col=colors[ord])
arrows(xx, y-ysdord/2, xx, y+(ysdord/2), length=0.05, angle=90, code=3, lwd=.2)

# without RNR2
xx=barplot(y[1:10],ylim=c(0,18), names=genesord[1:10], las=2, cex.names=.7, ylab='Ribo unbound half-life (min)', col=colors[ord])

arrows(xx, y[1:10]-ysdord[1:10]/2, xx, y[1:10]+(ysdord[1:10]/2), length=0.05, angle=90, code=3, lwd=.2)


dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlot_RiboAssocComparison_matrix.R')

