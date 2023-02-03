library(data.table)
library('scales')


path = '/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/TE/'

fns = c('Fibro', 'HEK', 'Myocyte', 'Myoblast')

# Read in data
DT1 <- data.table(read.table(paste0(path, fns[1], '_TE.txt'), header=TRUE))
DT2 <- data.table(read.table(paste0(path, fns[2], '_TE.txt'), header=TRUE))
DT3 <- data.table(read.table(paste0(path, fns[3], '_TE.txt'), header=TRUE))
DT4 <- data.table(read.table(paste0(path, fns[4], '_TE.txt'), header=TRUE))

# Combine DTs
DT <- cbind(DT1, DT2[,2:3], DT3[,2:3], DT4[,2:3])
# Remove ND6
DT <- DT[GeneName !=  'MT-ND6']

# Rename columns
sampnames <- c('Fibroblast', 'HEK293T_1', 'HEK293T_2','Myocyte_1','Myocyte_2', 'Myoblast_1', 'Myoblast_2')
colnames(DT) <- c('GeneName',sampnames)


# Set colors
col1='darkgoldenrod1' # 
col2='dodgerblue' # CI
col3='deeppink3' # CIV
col4='grey60' # CV
col5='aquamarine2' # CIII
colors = c(rep(col2,6), rep(col5,1), rep(col3,3), rep(col4,2))


ptsize = 1.5

pdf(paste0(path, 'TECorrelationMatrix.pdf'), 
     width=5.5, # width=4.5, width=2, width=2.5
     height=5.5, #4
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


# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlot_TEComparison_matrix.R')

