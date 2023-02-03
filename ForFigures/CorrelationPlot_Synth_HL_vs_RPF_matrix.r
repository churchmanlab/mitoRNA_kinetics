library(data.table)
library('scales')
library('rlist')
library(gclus)

measure='HL' # HL ProtSynth
path = '/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/CombinedData/'
filename = 'Gene_expression_model2.txt' # Gene_expression_model2.txt  All_abundances_2.txt

# Read in data
DTin <- data.table(read.table(paste0(path, filename), header=TRUE, skip=2, sep = '\t'))

cells=c('HelaModel', 'HelaRep2','HEKwt','HEKLRPPRC','K562') 
xCol='RPF' # HelaModel_RPF HelaRep2_RPF HEKwt_RPF HEKLRPPRC_RPF K562_RPF
yCol=measure 

xcols=paste0(cells, '_', xCol)
ycols=paste0(cells, '_', yCol)
cols=c('Gene',xcols, ycols)
DT = DTin[,..cols]


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
ord_index = append(ord_index, which(grepl(paste0('^',txpt,'$'), DT$Gene)))
}
setorder(DT[, .r := order(ord_index)], .r)[, .r := NULL]

ptsize = 1.5

pdf(paste0(path, measure, '_vs_RPF_CorrelationMatrix.pdf'), 
     width=8, # width=4.5, width=2, width=2.5
     height=8, #4
     pointsize=11)

# Correlation in absolute terms for coloring plots
corr <- abs(cor(DT[,2:length(DT)])) 
bgcol <- dmat.color(corr)
ord <- c(1:10) # order.single(corr) c(1:10)

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
  
#   cpairs(data,                    # Data frame of variables
#        order,                   # Order of the variables
#        panel.colors = colors,   # Matrix of panel colors
#        border.color = "grey70", # Borders color
#        gap = 0.45,              # Distance between subplots
#        main = "Ordered variables colored by correlation", # Main title
#        show.points = TRUE,      # If FALSE, removes all the points
#        pch = 21,                # pch symbol
#        bg = rainbow(3)[iris$Species]) # Colors by group

# Create plots
pairs(DT[,2:length(DT)],lower.panel = DataPanel, upper.panel = NULL, gap=.4, labels=colnames(DT[,2:length(DT)]))


cpairs(DT[,2:length(DT)],ord, panel.colors <- bgcol, gap=.4,labels=colnames(DT[,2:length(DT)]))

dev.off()
# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/CorrelationPlot_Synth_HL_vs_RPF_matrix.R')

