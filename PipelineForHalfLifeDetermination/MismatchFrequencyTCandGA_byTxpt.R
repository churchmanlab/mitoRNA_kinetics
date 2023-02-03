####!/n/app/R/4.0.1/bin/Rscript

###
# First make shortened file with run_AWKforTCperT_GAperG_fragments
 
# Top strand: 99=R1+, 147=R2- (primary), 355=R1+, 403=R2- (not primary)

### USE:
# sbatch -p short -t 0-01:00 --mem=50G --wrap="Rscript ../Scripts/MismatchFrequencyTCandGA_byTxpt.R"

library(data.table)
library(stringr)
library(purrr)
library(rlist)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

Exp = args[1]
Libs = unlist(strsplit(args[2],","))

MapMethod <- args[3]
reads <- args[4] # 'MTall' All MTnorRNA

path = paste0(getwd(), '/') 

samples=paste0(Exp, '_', Libs)
set = args[5] #'norRNA' # norRNA LowNew
task = args[6] # plotOnly All

NumSamps = length(samples)

shortnames = samples


TcountsEachRegion <- c()
GcountsEachRegion <- c()
TCs0EachRegion <- c()
GAs0EachRegion <- c()
TCs1EachRegion <- c()
GAs1EachRegion <- c()
TCs2EachRegion <- c()
GAs2EachRegion <- c()


TxptsAll = c('RNR2', 'ND1', 'ND2', 'CO1','CO2','ATP8_6','CO3','ND3','ND4L_4','ND5','CYB','7S','antiCYB','ND6', 'antiND5','antiND4L_4','antiND3', 'antiATP8_6_CO3', 'antiCO2','antiCO1', 'antiND2', 'antiND1')
Txpts = c('ND1', 'ND2', 'CO1','CO2','ATP8_6','CO3','ND3','ND4L_4','ND5','CYB','antiCYB','ND6', 'antiND5','antiND4L_4','antiND3', 'antiATP8_6_CO3', 'antiCO2','antiCO1', 'antiND2', 'antiND1', '7S')

Txpts2 = c('RNR2', 'CO1', 'CO2', 'ATP8_6', 'CO3', 'ND4L_4', 'antiCYB', 'ND6', 'antiND5', 'antiND4L_4', 'antiATP8_6_CO3', 'antiND1')

if (set == 'all') {
TOI = TxptsAll
}
if (set == 'norRNA') {
TOI = Txpts
}
if (set == 'LowNew') {
TOI = Txpts2
}


if (task != 'plotOnly') {

for (region in TOI) {

# Get data, make into DTs
for (i in c(1:NumSamps)) {
# get reads
assign(paste0('mDT',i), data.table(read.table(paste0(path, samples[i],'_',MapMethod,'/frag_MMfrequency_',reads,'.txt'), header=TRUE, sep='\t', quote='',stringsAsFactors = FALSE)))

# chr start end strand seq frag_length tot_mismatches TC_mismatches GA_mismatches
# Keep only reads in region


RNR2 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$end < 3229
ND1 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 3303 & get(paste0('mDT',i))$end < 4262
antiND1 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 3227 & get(paste0('mDT',i))$end < 4328
ND2 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 4468 & get(paste0('mDT',i))$end < 5511
antiND2 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 4399 & get(paste0('mDT',i))$end < 5586
CO1 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 5899 & get(paste0('mDT',i))$end < 7445
antiCO1 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 5897 & get(paste0('mDT',i))$end < 7445
CO2 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 7584 & get(paste0('mDT',i))$end < 8291
antiCO2 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 7513 & get(paste0('mDT',i))$end < 8276
ATP8_6 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 8363 & get(paste0('mDT',i))$end < 9207
CO3 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 9205 & get(paste0('mDT',i))$end < 9990
antiATP8_6_CO3 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 8283 & get(paste0('mDT',i))$end < 9990
ND3 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 10057 & get(paste0('mDT',i))$end < 10404
antiND3 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 10057 & get(paste0('mDT',i))$end < 10456
ND4L_4 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 10468 & get(paste0('mDT',i))$end < 12137
antiND4L_4 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 10468 & get(paste0('mDT',i))$end < 12137
ND5 = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 14741
antiND5 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 14111
ND6 = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14110 & get(paste0('mDT',i))$end < 14673
CYB = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 15887
antiCYB = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14745 & get(paste0('mDT',i))$end < 15951
antiCYB = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14745 & get(paste0('mDT',i))$end < 15951
MT7S = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 209 & get(paste0('mDT',i))$end < 409



if (region == 'all') {assign(paste0('mDT',i), get(paste0('mDT',i)))}
if (region == 'RNR2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[RNR2])
}
if (region == 'ND1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 3303 & get(paste0('mDT',i))$end < 4262])
}
if (region == 'antiND1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 3227 & get(paste0('mDT',i))$end < 4328])
}
if (region == 'ND2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 4468 & get(paste0('mDT',i))$end < 5511])
}
if (region == 'antiND2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 4399 & get(paste0('mDT',i))$end < 5586])
}
if (region == 'CO1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 5899 & get(paste0('mDT',i))$end < 7445])
}
if (region == 'antiCO1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 5897 & get(paste0('mDT',i))$end < 7445]) #  & get(paste0('mDT',i))$end !< 7825]
}
if (region == 'CO2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 7584 & get(paste0('mDT',i))$end < 8291])
}
if (region == 'antiCO2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 7513 & get(paste0('mDT',i))$end < 8276])
}
if (region == 'ATP8_6') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 8363 & get(paste0('mDT',i))$end < 9207])
}
if (region == 'CO3') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 9205 & get(paste0('mDT',i))$end < 9990])
}
if (region == 'antiATP8_6_CO3') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 8283 & get(paste0('mDT',i))$end < 9990])
}
if (region == 'ND3') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 10057 & get(paste0('mDT',i))$end < 10404])
}
if (region == 'antiND3') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 10057 & get(paste0('mDT',i))$end < 10456])
}
if (region == 'ND4L_4') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 10468 & get(paste0('mDT',i))$end < 12137])
}
if (region == 'antiND4L_4') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 10468 & get(paste0('mDT',i))$end < 12137])
}
if (region == 'ND5') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 14741])
}
if (region == 'antiND5') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 14111])
}
if (region == 'ND6') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14110 & get(paste0('mDT',i))$end < 14673])
}
if (region == 'CYB') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 12335 & get(paste0('mDT',i))$end < 15887])
}
if (region == 'antiCYB') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14745 & get(paste0('mDT',i))$end < 15951])
}


# Make new columns with number of Ts, Gs and the conversion rates
get(paste0('mDT',i))[, Tcount := ifelse(strand == '+', as.numeric(lapply(str_count(seq, pattern='T'),'[[',1)) + TC_mismatches, as.numeric(lapply(str_count(seq, pattern='A'),'[[',1)) + TC_mismatches)]

get(paste0('mDT',i))[, Gcount := ifelse(strand == '+', as.numeric(lapply(str_count(seq, pattern='G'),'[[',1)) + GA_mismatches, as.numeric(lapply(str_count(seq, pattern='C'),'[[',1)) + GA_mismatches)]
# get(paste0('nDT',i))[, ConvRate := round(TC_mismatches/Tcount, digits = 2)]
get(paste0('mDT',i))[, ConvRateTC := round(TC_mismatches/Tcount, digits = 2)]
get(paste0('mDT',i))[, ConvRateGA := round(GA_mismatches/Gcount, digits = 2)]


# Get average and std dev of T count and G count

assign(paste0('meanTcount',i), format(mean(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))
assign(paste0('Tstddev',i), format(sd(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))

assign(paste0('meanGcount',i), format(mean(get(paste0('mDT',i))$Gcount, na.rm=TRUE), digits=3))
assign(paste0('Gstddev',i), format(sd(get(paste0('mDT',i))$Gcount, na.rm=TRUE), digits=3))

# Get number of reads with x TC conv (use 30 as arbitrary max)
for (j in c(0:30)){
assign(paste0(paste0('TC',j,'_'),i),nrow(get(paste0('mDT',i))[get(paste0('mDT',i))$TC_mismatches == j]))
}

for (j in c(0:30)){
assign(paste0(paste0('GA',j,'_'),i),nrow(get(paste0('mDT',i))[get(paste0('mDT',i))$GA_mismatches == j]))
}

# Get total number of reads 
assign(paste0('tot',i),nrow(get(paste0('mDT',i))))

# Get fragment lengths and T counts and G counts
assign(paste0('lengths', i),nchar(get(paste0('mDT',i))$seq))
assign(paste0('TcountsDist', i),get(paste0('mDT',i))$Tcount)
assign(paste0('GcountsDist', i),get(paste0('mDT',i))$Gcount)
}




# Plot average # Ts and Gs per read with std dev and num with TC and GA MM

# Get max number of TC conversions in data
TC_mm = c()
GA_mm = c()
for (i in c(1:NumSamps)) {
TC_mm = c(TC_mm, get(paste0('mDT',i))$TC_mismatches)
GA_mm = c(GA_mm, get(paste0('mDT',i))$GA_mismatches)
}
maxTC = max(TC_mm)
maxGA = max(GA_mm)


Tcountslist=c()
Tsdevlist=c()
Gcountslist=c()
Gsdevlist=c()
totslist=c()

for (i in c(1:NumSamps)){
Tcountslist=c(Tcountslist, get(paste0('meanTcount',i)))
Tsdevlist=c(Tsdevlist, get(paste0('Tstddev',i)))
Gcountslist=c(Gcountslist, get(paste0('meanGcount',i)))
Gsdevlist=c(Gsdevlist, get(paste0('Gstddev',i)))
totslist=c(totslist, get(paste0('tot',i)))
}

Tcounts <- as.numeric(Tcountslist) 
Tsdevs <- as.numeric(Tsdevlist) 
Gcounts <- as.numeric(Gcountslist) 
Gsdevs <- as.numeric(Gsdevlist) 
tots <-  as.numeric(totslist)
minReadCounts = min(totslist[1:(length(totslist))])

ylimitsCounts = c(0,55)

TcountsEachRegion <- rbind(TcountsEachRegion, Tcounts) # mean Ts per read in the region
GcountsEachRegion <- rbind(GcountsEachRegion, Gcounts)

for (j in c(0:maxTC)){
assign(paste0('TClist',j),c())
assign(paste0('GAlist',j),c())
}

# Number of reads with j TC or GA mismatches in each sample (i) in the region
for (j in c(0:maxTC)){
for (i in c(1:NumSamps)){
assign(paste0('TClist',j),c(get(paste0('TClist',j)), get(paste0(paste0('TC',j,'_'),i)))) 
assign(paste0('GAlist',j),c(get(paste0('GAlist',j)), get(paste0(paste0('GA',j,'_'),i))))
}
}

# For frequency of each count
for (j in c(0:maxTC)){
assign(paste0('TC',j,'s'),as.numeric(get(paste0('TClist',j))))
assign(paste0('GA',j,'s'),as.numeric(get(paste0('GAlist',j))))
}

Tvec = c()
Gvec = c()
for (j in c(0:maxTC)){
Tvec = c(Tvec,get(paste0('TC',j,'s')))
Gvec = c(Gvec,get(paste0('GA',j,'s')))
}

# Make matrix where rows are number of TC/GA mismatches (1,2,3,..max), columns are samples, entries are number of reads with that many mismatches
TCs <- matrix(Tvec, ncol=NumSamps, nrow=maxTC+1, byrow=TRUE) 
GAs <- matrix(Gvec, ncol=NumSamps, nrow=maxTC+1, byrow=TRUE)
readcountMat = matrix(tots, ncol=NumSamps, nrow=maxTC+1, byrow=TRUE)

# For frequency of each count, normalized to per 1000 reads
# Make matrix with readcounts to divide TCs matrix
NormTCs = round(TCs/readcountMat*1000, 0)
NormGAs = round(GAs/readcountMat*1000, 0)

TCs0 = NormTCs[1,]
GAs0 = NormGAs[1,]
TCs1 = NormTCs[2,]
GAs1 = NormGAs[2,]
TCs2 = NormTCs[3,]
GAs2 = NormGAs[3,]

TCs0EachRegion <- rbind(TCs0EachRegion, TCs0)
GAs0EachRegion <- rbind(GAs0EachRegion, GAs0)
TCs1EachRegion <- rbind(TCs1EachRegion, TCs1)
GAs1EachRegion <- rbind(GAs1EachRegion, GAs1)
TCs2EachRegion <- rbind(TCs2EachRegion, TCs2)
GAs2EachRegion <- rbind(GAs2EachRegion, GAs2)
}


write.table(TCs1EachRegion, file=paste0(path, 'MMfrequency/TCs1EachRegion.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(GAs1EachRegion, file=paste0(path, 'MMfrequency/GAs1EachRegion.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(TCs2EachRegion, file=paste0(path, 'MMfrequency/TCs2EachRegion.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(GAs2EachRegion, file=paste0(path, 'MMfrequency/GAs2EachRegion.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)







pdf(paste0(path, 'MMfrequency/',Exp,'_', MapMethod, '_frag_AveTsGs_NumReadsWithMM_eachTxpt_', set,'.pdf'), width = 16, height = 18)

# par(mfrow=c(3,4), cex.lab = 1)
par(mfrow=c(6,4), cex.lab = 1)

lineweight = 2.5


########################################
# Plot T and G counts
for (i in c(1:NumSamps)) {
plot(TcountsEachRegion[,i], type='l', col='red', main = paste0(shortnames[i], ' T/G counts'), ylab = 'Mean counts per read', ylim = c(0, 45), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
lines(GcountsEachRegion[,i], type='l', col='orange',lwd=lineweight)
legend('bottomright', legend=c('T', 'G'), bty="n", text.col =  c('red', 'orange'))
}

########################################
# Plot 1 TC and GA per read across txpts
for (i in c(1:NumSamps)) {
plot(TCs1EachRegion[,i], type='l', col='gold1', main = shortnames[i], ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 300), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
lines(GAs1EachRegion[,i], type='l', col='dodgerblue', lwd=lineweight)
legend('topright', legend=c('TC (1/read)', 'GA (1/read)'), bty="n", text.col =  c('gold1', 'dodgerblue'))
}

########################################

if (Exp == 'TL6') {
cols = c('grey60', 'gold1', 'dodgerblue', 'forestgreen')
bg=1
sU4_sl=2
sG6_sl200=3
dl_200=4
}
if (Exp == 'TL7') {
cols = c('grey60', 'gold1', 'dodgerblue', 'darkblue', 'yellowgreen', 'forestgreen')
bg=1
sU4_sl=2
sG6_sl200=3
sG6_sl1000=4
dl_200=5
dl_1000=6
}

# Plot 1 TC per read across txpts, for ALL samples
plot(TCs1EachRegion[,1], type='l', col=cols[1], main = '1 TC/read', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 350), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
for (i in c(2:NumSamps)) {
lines(TCs1EachRegion[,i], type='l', col=cols[i], lwd=lineweight)
}
legend('topright', legend=shortnames, bty="n", text.col =  cols)

# Plot 1 GA per read across txpts, for ALL samples
plot(GAs1EachRegion[,1], type='l', col=cols[1], main = '1 GA/read', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 60), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
for (i in c(2:NumSamps)) {
lines(GAs1EachRegion[,i], type='l', col=cols[i], lwd=lineweight)
}
legend('topleft', legend=shortnames, bty="n", text.col =  cols)

# Plot 2 TC per read across txpts, for ALL samples
plot(TCs2EachRegion[,1], type='l', col=cols[i], main = '2 TC/read', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 120), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
for (i in c(2:NumSamps)) {
lines(TCs2EachRegion[,i], type='l', col=cols[i], lwd=lineweight)
}
legend('topright', legend=shortnames, bty="n", text.col =  cols)

# Plot 2 GA per read across txpts, for ALL samples
plot(GAs2EachRegion[,1], type='l', col=cols[i], main = '2 GA/read', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 10), lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2)
for (i in c(2:NumSamps)) {
lines(GAs2EachRegion[,i], type='l', col=cols[i], lwd=lineweight)
}
legend('topleft', legend=shortnames, bty="n", text.col =  cols)


########################################
#### Have to manually modify this part
########################################

# Norm to single label, NO background subtraction
# ylimits = c(0,2.5)
# 
# 
# plot(TCs1EachRegion[,5]/TCs1EachRegion[,2], type='l', col=cols[5], main = '1 TC/read\nNorm to single label', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# lines(TCs1EachRegion[,6]/TCs1EachRegion[,2], type='l', col=cols[6], lwd=lineweight)
# legend('topright', legend=c('180m/90m 4sU (200 6sG)', '180m/90m 4sU (1000 6sG)'), bty="n", text.col =  cols[5:6])
# abline(h=1, lty = 2)
# 
# plot(GAs1EachRegion[,5]/GAs1EachRegion[,3], type='l', col=cols[5], main = '1 GA/read\nNorm to single label', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim=ylimits,lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# lines(GAs1EachRegion[,6]/GAs1EachRegion[,4], type='l', col=cols[6], lwd=lineweight)
# legend('topleft', legend=c('90m dl/90m 6sG 200','90m dl/90m 6sG 1000'), bty="n", text.col =  cols[5:6])
# abline(h=1, lty = 2)
# 
# plot(TCs2EachRegion[,5]/TCs2EachRegion[,2], type='l', col=cols[5], main = '2 TC/read\nNorm to single label', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)',ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# lines(TCs2EachRegion[,6]/TCs2EachRegion[,2], type='l', col=cols[6], lwd=lineweight)
# legend('topright', legend=c('180m/90m 4sU (200 6sG)', '180m/90m 4sU (1000 6sG)'), bty="n", text.col =  cols[5:6])
# abline(h=1, lty = 2)
# 
# plot(GAs2EachRegion[,5]/GAs2EachRegion[,3], type='l', col=cols[5], main = '2 GA/read\nNorm to single label', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)',ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# lines(GAs2EachRegion[,6]/GAs2EachRegion[,4], type='l', col=cols[6], lwd=lineweight)
# legend('topleft', legend=c('90m dl/90m 6sG 200','90m dl/90m 6sG 1000'), bty="n", text.col =  cols[5:6])
# abline(h=1, lty = 2)
# 
#######################################
### BACKGROUND SUBTRACTED
########################################

# Plot 1 TC per read across txpts, for ALL samples, background subtracted
# plot(TCs1EachRegion[,1]-TCs1EachRegion[,1], type='l', col=cols[1], main = '1 TC/read, background subtracted', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 350), lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# for (i in c(2:NumSamps)) {
# lines(TCs1EachRegion[,i]-TCs1EachRegion[,1], type='l', col=cols[i], lwd=lineweight)
# }
# legend('topright', legend=shortnames, bty="n", text.col = cols)
# 
# # Plot 1 GA per read across txpts, for ALL samples
# plot(GAs1EachRegion[,1]-GAs1EachRegion[,1], type='l', col=cols[1], main = '1 GA/read, background subtracted', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 60), lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# for (i in c(2:NumSamps)) {
# lines(GAs1EachRegion[,i]-GAs1EachRegion[,1], type='l', col=cols[i], lwd=lineweight)
# }
# legend('topleft', legend=shortnames, bty="n", text.col = cols)
# 
# # Plot 2 TC per read across txpts, for ALL samples
# plot(TCs2EachRegion[,1]-TCs2EachRegion[,1], type='l', col=cols[1], main = '2 TC/read, background subtracted', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 120), lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# for (i in c(2:NumSamps)) {
# lines(TCs2EachRegion[,i]-TCs2EachRegion[,1], type='l', col=cols[i], lwd=lineweight)
# }
# legend('topright', legend=shortnames, bty="n", text.col = cols)
# 
# # Plot 2 GA per read across txpts, for ALL samples
# plot(GAs2EachRegion[,1]-GAs2EachRegion[,1], type='l', col=cols[1], main = '2 GA/read, background subtracted', ylab = 'Normalized frequency (counts/1k reads mapped to txpt)', ylim = c(0, 10), lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# for (i in c(2:NumSamps)) {
# lines(GAs2EachRegion[,i]-GAs2EachRegion[,1], type='l', col=cols[i], lwd=lineweight)
# }
# legend('topleft', legend=shortnames, bty="n", text.col = cols)
# 

########################################

ylimits = c(-3,3)
ptsize = 1.5

## For removing ND3
# TOI = TOI[-7] 
# TCs1EachRegion = TCs1EachRegion[-7,]
# GAs1EachRegion = GAs1EachRegion[-7,]
# TCs2EachRegion = TCs2EachRegion[-7,]
# GAs2EachRegion = GAs2EachRegion[-7,]
Txpts = c('RNR2', 'ND1', 'ND2', 'CO1','CO2','ATP8_6','CO3','ND3','ND4L_4','ND5','CYB','antiCYB','ND6', 'antiND5','antiND4L_4','antiND3', 'antiATP8_6_CO3', 'antiCO2','antiCO1', 'antiND2', 'antiND1', '7S')



if (set == 'all') {
heavy=1:11
light=12:length(TOI)}
if (set == 'norRNA') {
heavy=1:10
light=11:length(TOI)}
if (set == 'LowNew') {
heavy=1:6
light=7:length(TOI)}

nheavy = max(heavy)
ntot = max(light)

# Plot 1TC/read background subtracted,  norm to single label (200 6sG)
Norm = (TCs1EachRegion[,dl_200]-TCs1EachRegion[,bg])/(TCs1EachRegion[,sU4_sl]-TCs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy])
HsAve = mean(Norm[heavy])

HSline = HsAve
HslineLegend = 'Heavy strand (mean)'
LSline = LsAve
LslineLegend = 'Light strand (mean)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='orange', main = '1 TC/read, background subtracted\nNorm to single label (200 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('180m/90m 4sU (200 6sG)', HslineLegend, LslineLegend), bty="n", text.col =  c('orange', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

if (Exp == 'TL7') {
# Plot 1TC/read background subtracted,  norm to single label (1000 6sG)
Norm = (TCs1EachRegion[,dl_1000]-TCs1EachRegion[,bg])/(TCs1EachRegion[,sU4_sl]-TCs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy])
HsAve = mean(Norm[heavy])

HSline = HsAve
HslineLegend = 'Heavy strand (mean)'
LSline = LsAve
LslineLegend = 'Light strand (mean)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='orange', main = '1 TC/read, background subtracted\nNorm to single label (1000 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('180m/90m 4sU (1000 6sG)', HslineLegend, LslineLegend), bty="n", text.col =  c('orange', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

}


# Plot 1GA/read background subtracted,  norm to single label (200 6sG)
Norm = (GAs1EachRegion[,dl_200]-GAs1EachRegion[,bg])/(GAs1EachRegion[,sG6_sl200]-GAs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)
HsAve = mean(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)

HSline = HsAve
HslineLegend = 'Heavy strand (mean)'
LSline = LsAve
LslineLegend = 'Light strand (mean)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='blue', main = '1 GA/read, background subtracted\nNorm to single label (200 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits,lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('90m dl/90m 6sG (200)', HslineLegend, LslineLegend), bty="n", text.col =  c('blue', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

if (Exp == 'TL7') {
# Plot 1GA/read background subtracted,  norm to single label (1000 6sG)
Norm = (GAs1EachRegion[,6]-GAs1EachRegion[,bg])/(GAs1EachRegion[,4]-GAs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)
HsAve = mean(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)

HSline = HsAve
HslineLegend = 'Heavy strand (mean)'
LSline = LsAve
LslineLegend = 'Light strand (mean)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='blue', main = '1 GA/read, background subtracted\nNorm to single label (1000 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits,lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('90m dl/90m 6sG (1000)', HslineLegend, LslineLegend), bty="n", text.col =  c('blue', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)
}

dev.off()


} # End of region to skip if plotOnly



if (task == 'plotOnly') {


TCs1EachRegion <- as.matrix(read.table(paste0(path, 'MMfrequency/TCs1EachRegion.txt'), header=FALSE, sep='\t', quote='',stringsAsFactors = FALSE))
GAs1EachRegion <- as.matrix(read.table(paste0(path, 'MMfrequency/GAs1EachRegion.txt'), header=FALSE, sep='\t', quote='',stringsAsFactors = FALSE))
TCs2EachRegion <- as.matrix(read.table(paste0(path, 'MMfrequency/TCs2EachRegion.txt'), header=FALSE, sep='\t', quote='',stringsAsFactors = FALSE))
GAs2EachRegion <- as.matrix(read.table(paste0(path, 'MMfrequency/GAs2EachRegion.txt'), header=FALSE, sep='\t', quote='',stringsAsFactors = FALSE))

if (Exp == 'TL6') {
cols = c('grey60', 'gold1', 'dodgerblue', 'forestgreen')
bg=1
sU4_sl=2
sG6_sl200=3
dl_200=4
}
if (Exp == 'TL7') {
cols = c('grey60', 'gold1', 'dodgerblue', 'darkblue', 'yellowgreen', 'forestgreen')
bg=1
sU4_sl=2
sG6_sl200=3
sG6_sl1000=4
dl_200=5
dl_1000=6
}

if (set == 'all') {
heavy=1:11
light=12:length(TOI)}
if (set == 'norRNA') {
heavy=1:10
light=11:length(TOI)}
if (set == 'LowNew') {
heavy=1:6
light=7:length(TOI)}

nheavy = max(heavy)
ntot = max(light)

ylimits=c(-3,3)
ptsize = 1.5
lineweight = 2.5
}



# Plot exact same thing again, as its own pdf

pdf(paste0(path, 'MMfrequency/',Exp,'_', MapMethod, '_DoublabelNormToSingle_eachTxpt_', set,'.pdf'), width = 5, height = 4.5)

# Plot 1TC/read background subtracted,  norm to single label (200 6sG)
Norm = (TCs1EachRegion[,dl_200]-TCs1EachRegion[,bg])/(TCs1EachRegion[,sU4_sl]-TCs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy])
HsAve = mean(Norm[heavy])

HSline = HsMed
HslineLegend = 'Heavy strand (median)'
LSline = LsMed
LslineLegend = 'Light strand (median)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='orange', main = '1 TC/read, background subtracted\nNorm to single label (200 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('180m/90m 4sU (200 6sG)', HslineLegend, LslineLegend), bty="n", text.col =  c('orange', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

if (Exp == 'TL7') {
# Plot 1TC/read background subtracted,  norm to single label (1000 6sG)
Norm = (TCs1EachRegion[,dl_1000]-TCs1EachRegion[,bg])/(TCs1EachRegion[,sU4_sl]-TCs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy])
HsAve = mean(Norm[heavy])

HSline = HsMed
HslineLegend = 'Heavy strand (median)'
LSline = LsMed
LslineLegend = 'Light strand (median)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='orange', main = '1 TC/read, background subtracted\nNorm to single label (1000 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('180m/90m 4sU (1000 6sG)', HslineLegend, LslineLegend), bty="n", text.col =  c('orange', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

}

# Plot 1GA/read background subtracted,  norm to single label (200 6sG)
Norm = (GAs1EachRegion[,dl_200]-GAs1EachRegion[,bg])/(GAs1EachRegion[,sG6_sl200]-GAs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)
HsAve = mean(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)

HSline = HsMed
HslineLegend = 'Heavy strand (median)'
LSline = LsMed
LslineLegend = 'Light strand (median)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='blue', main = '1 GA/read, background subtracted\nNorm to single label (200 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits,lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('90m dl/90m 6sG (200)', HslineLegend, LslineLegend), bty="n", text.col =  c('blue', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

if (Exp == 'TL7') {
# Plot 1GA/read background subtracted,  norm to single label (1000 6sG)
Norm = (GAs1EachRegion[,6]-GAs1EachRegion[,bg])/(GAs1EachRegion[,4]-GAs1EachRegion[,bg])

LsMed = median(Norm[light])
LsAve = mean(Norm[light])
HsMed = median(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)
HsAve = mean(Norm[heavy][is.finite(Norm[heavy])], na.rm=TRUE)

HSline = HsMed
HslineLegend = 'Heavy strand (median)'
LSline = LsMed
LslineLegend = 'Light strand (median)'

plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='blue', main = '1 GA/read, background subtracted\nNorm to single label (1000 6sG)', ylab = 'Normalized frequency (log2)', ylim=ylimits,lwd=lineweight, xaxt = 'n', xlab = ' ')
axis(1, at=1:length(TOI), labels=TOI, las=2, cex=.8)
legend('topright', legend=c('90m dl/90m 6sG (1000)', HslineLegend, LslineLegend), bty="n", text.col =  c('blue', 'black', 'red'))
abline(h=0, lty = 1, lwd = .5,col = 'red')
abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
polygon(c(0,0,nheavy+.5,nheavy+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
polygon(c(nheavy+.5,nheavy+.5,ntot+.5,ntot+.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)
}



# Left from TL6, for plotting 2 TC or GA/read (not enough of these for anything to show up)
# Norm = (TCs2EachRegion[,4]-TCs2EachRegion[,1])/(TCs2EachRegion[,2]-TCs2EachRegion[,1])
# 
# LsMed = median(Norm[light])
# LsAve = mean(Norm[light])
# HsMed = median(Norm[heavy])
# HsAve = mean(Norm[heavy])
# 
# HSline = HsAve
# HslineLegend = 'Heavy strand (mean)'
# LSline = LsAve
# LslineLegend = 'Light strand (mean)'
# 
# plot(log(Norm/LSline,2), type='p',pch = 16, cex=ptsize, col='orange', main = '2 TC/read, background subtracted\nNorm to single label', ylab = 'Normalized frequency (log2)',ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# legend('topright', legend=c('180m/90m 4sU', HslineLegend, LslineLegend), bty="n", text.col =  c('orange', 'black', 'red'))
# abline(h=0, lty = 1, lwd = .5,col = 'red')
# abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
# polygon(c(0,0,9.5,9.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
# polygon(c(9.5,9.5,20,20), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)
# 
# Norm = (GAs2EachRegion[,4]-GAs2EachRegion[,1])/(GAs2EachRegion[,3]-GAs2EachRegion[,1])
# 
# LsMed = median(Norm[light])
# LsAve = mean(Norm[light])
# HsMed = median(Norm[light])
# HsAve = mean(Norm[light])
# 
# HSline = HsAve
# HslineLegend = 'Heavy strand (mean)'
# LSline = LsAve
# LslineLegend = 'Light strand (mean)'
# 
# plot(log(Norm/LSline,2),  type='p',pch = 16, cex=ptsize, col='blue', main = '2 GA/read, background subtracted\nNorm to single label', ylab = 'Normalized frequency (log2)',ylim=ylimits, lwd=lineweight, xaxt = 'n', xlab = ' ')
# axis(1, at=1:length(TOI), labels=TOI, las=2)
# legend('topleft', legend=c('180m/90m 6sG'), bty="n", text.col =  c('blue'))
# abline(h=0, lty = 1, lwd = .5,col = 'red')
# abline(h=log(HSline/LSline,2), lty = 2, col = 'black')
# polygon(c(0,0,9.5,9.5), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('black', .1), border = NA)
# polygon(c(9.5,9.5,20,20), c(min(ylimits),max(ylimits),max(ylimits),min(ylimits)), col = alpha('red', .1), border = NA)

#######################################
dev.off()

# source('/n/groups/churchman/mc348/TimelapseSeq/Scripts/MismatchFrequencyTCandGA_byTxpt.R')
