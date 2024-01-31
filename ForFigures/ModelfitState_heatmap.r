library(data.table)
library('scales')


Set = 'tot' # tot mito_enrich


DT <- data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/HalfLife/TL4_MT_t5MTMMinformed6_modeAll_PcMTnorRNA_FracNew_halflives_corr_1592min_modelprobs.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))
# Processing_kinetics_w_mito-enriched.txt has mito-enriched sample pooled with the total

if (Set == 'tot') {
samp = 'Half.life.unprocessed..min.'
} else if (Set == 'mito_enrich') {
samp = 'Half.life.unprocessed..min._mito.enriched'
}

pdf(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Processing/HLvsConNoncon_', Set,'.pdf'), width = 3, height = 4)

# Make a new column with each category as 1,2,3,4 to plot them in the right place
PK[, x := ifelse(Canonical == 'Yes' & Transcript_end=='5 prime', 1, ifelse(Canonical == 'Yes' & Transcript_end=='3 prime', 2, ifelse(Canonical == 'No' & Transcript_end=='5 prime', 3, 4)))]

x=PK$x
y=c(PK[[samp]])

trans=.5
# Start plot without adding points
plot(x,y, type='n', xlim = c(0, 5), xaxt='n', xlab='', ylab="Half-life (min) of unprocessed", main=Set)
axis(1, at=c(1,2, 3,4), labels=c("Canonical 5'", "Canonical 3'","Non-can 5'", "Non-can 3'"), las=2)
# xx = boxplot(c(Con, NonCon), ylab = "Half-life (min)", names= c('Canonical', 'Non-canonical'), las=2, border='white')
# # Plot the data again into the same plot and customise the point shape, etc to your liking
points(jitter(x[PK$Transcript_end=='5 prime'],.5), y[PK$Transcript_end=='5 prime'], col=alpha('black', trans), pch=16, cex=1.5)
points(jitter(x[PK$Transcript_end=='3 prime'],.5), y[PK$Transcript_end=='3 prime'], col=alpha('red', trans), pch=16, cex=1.5)
# points(rep(2,nrow(NonCon)), NonCon[[1]], pch=16)
legend('topleft', c("5'", "3'"), text.col = c('black', 'red'), box.lty=0)


dev.off()

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/Processing_binaryScatter.r')

