library(data.table)
library('scales')


wid = 3.5 # width of plot 2.4
ht = 4 # height of plot
txt = .8 # axis text
las = 2



NucCol = "Average" # "Fraction_in_Ribo_ave" 
MitoCols = c('TL8', 'TL9','TL11')

DT1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/HalfLife/TL3_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_1592min_v2.txt'),header=TRUE, stringsAsFactors=FALSE, sep='\t'))

DT2 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_Nuc_t5MMinformed5_lenient_modeAll_FracNew_halflives_corr_1592min_v2.txt'),header=TRUE, stringsAsFactors=FALSE, sep='\t'))

# Get list of protein-coding genes
PCgenes <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/protein-coding_gene.txt'),header=TRUE, stringsAsFactors=TRUE, sep='\t', quote="", comment.char = ""))

symbol=as.vector(PCgenes$symbol)
alias=noquote(as.vector(PCgenes$alias_symbol))

namess=list()
for (i in c(1:length(alias))) {
names=strsplit(alias[i] , '["|]')[[1]]
namess[[length(namess) + 1]] <- names
}
aliases=unique(unlist(namess))

allnames = c(symbol, aliases)

# Get OXPHOS genes
OXPHOSgenes <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_list_newSymbols.txt'), header=FALSE, stringsAsFactors=FALSE))$V1


# Filter for only protein-coding genes
DT1ox <- DT1[DT1$Gene %in% OXPHOSgenes]
DT2ox <- DT2[DT2$Gene %in% OXPHOSgenes]

DT1oth <- DT1[DT1$Gene %in% allnames]
DT2oth <- DT2[DT2$Gene %in% allnames]

ox1 = DT1ox$Half.Life
ox2 = DT2ox$Half.Life
oth1 = DT1oth$Half.Life
oth2 = DT2oth$Half.Life


# Student's t-test
ttest_rep1 <- t.test(oth1, ox1)
ttest_rep2 <- t.test(oth2, ox2)

# Set which rows in xx boxplot table to use to remove outliers (1,5), or only use IQR (2,4)
row1=1
row2=5

ttest_rep1_noout <- t.test(oth1[oth1>xx$stats[row1,1] & oth1<xx$stats[row2,1]], ox1[ox1>xx$stats[row1,2] & ox1<xx$stats[row2,2]])

ttest_rep2_noout <- t.test(oth2[oth2>xx$stats[row1,3] & oth2<xx$stats[row2,3]], ox2[ox2>xx$stats[row1,4] & ox2<xx$stats[row2,4]])

p1 <- formatC(ttest_rep1_noout$p.value, format="e", digits=2)
p2 <- formatC(ttest_rep2_noout$p.value, format="e", digits=2)


# Make plot
pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/NucOther_OXPHOS_HL_Boxplot.pdf'), width = wid, height = ht)

xx=boxplot(oth1, ox1,oth2,ox2,  ylab = 'Half life (min)',  xlab = 'genes', names=c('Other rep1', 'OXPHOS rep1', 'Other rep2', 'OXPHOS rep2'), pch=16, outcex=.2, outcol=alpha('black', .2), cex.axis=txt,outline=FALSE,las=las, ylim=c(0, 1700))


text(1.5, max(xx$stats), label=p1)
text(3.5, max(xx$stats), label=p2)

dev.off()



# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/OXPHOS__allNuc_boxplots_HL.r')

