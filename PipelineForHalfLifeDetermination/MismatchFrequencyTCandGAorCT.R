####!/n/app/R/4.0.1/bin/Rscript

###
# First make shortened file with run_AWKforTCperT_GAperG_fragments
 
# Top strand: 99=R1+, 147=R2- (primary), 355=R1+, 403=R2- (not primary)

### USE:  
# Regions="ND1 antiND1 ND2 antiND2 CO1 antiCO1 CO2 antiCO2 ATP8_6 CO3 antiATP8_6_CO3 ND3 antiND3 ND4L_4 antiND4L_4 ND5 antiND5 ND6 CYB antiCYB" #  all/MTnorRNA/MTheavyNorRNA/MTlightNorRNA/LowNew/Bin1_2_3_4/RNR1"
# Regions="MTnorRNA"

library(data.table)
library(stringr)
library(purrr)
library(rlist)
args <- commandArgs(trailingOnly = TRUE)

Exp = args[1] # SStot1 NHC1
Libs = unlist(strsplit(args[2],","))
# Libs = c('0m','7m','15m','30m','45m','60m','90m','120m','240m')
# Libs = c('0m_tot','15m_tot','30m_tot','60m_tot', '0m_IP','15m_IP','30m_IP','60m_IP')
MapMethod <- args[3]
reads <- args[4] # 'MTall' All MTnorRNA
second <- args[5]
# region = 'MTnorRNA' # all MTnorRNA
region <- args[6]
stringentFilter <- args[7]
task <- args[8]

path = paste0(getwd(), '/') 
samples=paste0(Exp, '_', Libs)
seqmeth='paired'


NumSamps = length(samples)

shortnames = samples

# Get data, make into DTs
for (i in c(1:NumSamps)) {
# get reads
assign(paste0('mDT',i), data.table(read.table(paste0(path, samples[i],'_',MapMethod,'/frag_MMfrequency_',reads,'.txt'), header=TRUE, sep='\t', quote='',stringsAsFactors = FALSE)))

# chr start end strand seq frag_length tot_mismatches TC_mismatches GA_mismatches readName
# Keep only reads in region
x7S = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 209 & get(paste0('mDT',i))$end < 409
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
CYB = get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 14745 & get(paste0('mDT',i))$end < 15887
antiCYB = get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 14745 & get(paste0('mDT',i))$end < 15951



if (region == 'all') {assign(paste0('mDT',i), get(paste0('mDT',i)))}
if (region == 'MTnorRNA') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$start > 3229])
}
if (region == 'MTnoRNR1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$start > 1671])
}
if (region == 'MTsomerRNA') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$start > 2500]) # | get(paste0('mDT',i))$end < 1590])
}
if (region == 'RNR1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$end < 1590])
}
if (region == 'LowNew') {
assign(paste0('mDT',i), get(paste0('mDT',i))[RNR2 | antiATP8_6_CO3 | antiCO1 | antiCO2 | antiCYB | antiND1 | antiND4L_4 | antiND5 | ATP8_6 | CO1 | CO2 | CO3 | ND4L_4 | ND6 ])
}

if (region == 'MTheavyNorRNA') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$start > 3229])
}
if (region == 'MTlightNorRNA') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 3229])
}

if (region == 'RNR2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '+' & get(paste0('mDT',i))$end < 3229])
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
assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$strand == '-' & get(paste0('mDT',i))$start > 5897 & get(paste0('mDT',i))$end < 7445])
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
if (region == 'Bin1') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND2 | ND3 | CO3 ])
}
if (region == 'Bin2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND5 | ND6 | CYB ])
}
if (region == 'Bin3') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND4L_4 | ND1 ])
}
if (region == 'Bin4') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ATP8_6 | CO2 ])
}
if (region == 'Bin5') {
assign(paste0('mDT',i), get(paste0('mDT',i))[CO1])
}
if (region == 'allPC') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND1 | ND2 | ND3 | CO3 |ND5 | ND6 | CYB | ND4L_4 | ATP8_6 | CO2 | CO1])
}
if (region == 'Bin1_2') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND2 | ND3 | CO3 |ND5 | ND6 | CYB ])
}
if (region == 'Bin3_4') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND4L_4 | ND1 | ATP8_6 | CO2])
}
if (region == 'Bin1_2_3_4') {
assign(paste0('mDT',i), get(paste0('mDT',i))[ND2 | ND3 | CO3 |ND5 | ND6 | CYB | ND4L_4 | ND1 | ATP8_6 | CO2])
}


# remove fragments we determined to be misaligned
# `%notin%` <- Negate(`%in%`)
# assign(paste0('mDT',i), get(paste0('mDT',i))[get(paste0('mDT',i))$end %notin% c(2389, 2390, 3120,3121,3122,3123)])


if (task == 'readcountsOnly') {
# Get filtered read counts
lighttxpts=list(x7S,antiND1,antiND2,antiCO1,antiCO2,antiATP8_6_CO3,antiND3,antiND4L_4,antiND5,ND6,antiCYB)
lighttxptnames=c('7S','MT-antiND1', 'MT-antiND2','MT-antiCO1','MT-antiCO2','MT-antiATP8-6-CO3','antiND3','MT-antiND4L-4','MT-antiND5','MT-ND6','MT-antiCYB')
ltlengths = c(200,1100,1186,1547,762,1706,398,1668,1775,562,1196)

# Only "orange reads"
rcountsAll = c()
for (g in c(1:length(lighttxpts))) {
rcount = nrow(get(paste0('mDT',i))[lighttxpts[[g]]])
rcountsAll = c(rcountsAll, rcount)
}
rcountsTC = c()
for (g in c(1:length(lighttxpts))) {
rcount = nrow(get(paste0('mDT',i))[lighttxpts[[g]]][TC_mismatches > 0])
rcountsTC = c(rcountsTC, rcount)
}
# "Orange" reads and unlabelled
rcountsTCunlab = c()
for (g in c(1:length(lighttxpts))) {
rcount = nrow(get(paste0('mDT',i))[lighttxpts[[g]]][tot_mismatches - TC_mismatches < 1])
rcountsTCunlab = c(rcountsTCunlab, rcount)
}
RPK_All = rcountsAll/ltlengths*1000
RPK_TC = rcountsTC/ltlengths*1000
RPK_TCunlab = rcountsTCunlab/ltlengths*1000
# Normalize to antiND1
RPK_All_norm = RPK_All/RPK_All[2]
RPK_TC_norm = RPK_TC/RPK_TC[2]
RPK_TCunlab_norm = RPK_TCunlab/RPK_TCunlab[2]

dat = data.table(GeneName = lighttxptnames, TConlyRPK = RPK_TC, TCunlabRPK = RPK_TCunlab, TCallRPK = RPK_All, TConlyRPK_norm = RPK_TC_norm, TCunlabRPK_norm = RPK_TCunlab_norm, AllRPK_norm = RPK_All_norm)

write.table(dat, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_',samples[i],'_RPK.txt'), sep=("\t"), quote=FALSE, row.names=FALSE)
}


if (task != 'readcountsOnly') {
if (stringentFilter == 'yes') {
# remove fragments that have as many or more other mismatches as T>C
# allow up to n OTHER mismatches/100 nt
n=1
DT <- get(paste0('mDT',i))
filtDT <- DT[(DT$tot_mismatches - DT$TC_mismatches)/DT$frag_length*100 <= n]
# Write the read names
readNames=sapply(strsplit(filtDT$readName, '/'),"[[",1)
write.table(readNames, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_readNames_',samples[i],'_',n,'orLessOtherMM.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
assign(paste0('mDT',i), filtDT)
}



if (second == 'GA') {
ntF = 'G'
ntR = 'C'
ntconF = 'A'
ntconR = 'T'
mmcolumn = 'GA_mismatches'
} else if (second == 'CT') {
ntF = 'C'
ntR = 'G'
ntconF = 'T'
ntconR = 'A'
mmcolumn = 'CT_mismatches'
}


# Make new columns with number of Ts, Gs and the conversion rates
get(paste0('mDT',i))[, Tcount := ifelse(strand == '+', as.numeric(lapply(str_count(seq, pattern='T'),'[[',1)) + TC_mismatches, as.numeric(lapply(str_count(seq, pattern='A'),'[[',1)) + TC_mismatches)]

get(paste0('mDT',i))[, Ncount := ifelse(strand == '+', as.numeric(lapply(str_count(seq, pattern=ntF),'[[',1)) + get(mmcolumn), as.numeric(lapply(str_count(seq, pattern=ntR),'[[',1)) + get(mmcolumn))]

get(paste0('mDT',i))[, ConvRateTC := round(TC_mismatches/Tcount, digits = 2)]
get(paste0('mDT',i))[, ConvRateNN := round(get(mmcolumn)/Ncount, digits = 2)]



# Get average and std dev of T count and G count

assign(paste0('meanTcount',i), format(mean(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))
assign(paste0('Tstddev',i), format(sd(get(paste0('mDT',i))$Tcount, na.rm=TRUE), digits=3))

assign(paste0('meanNcount',i), format(mean(get(paste0('mDT',i))$Ncount, na.rm=TRUE), digits=3))
assign(paste0('Nstddev',i), format(sd(get(paste0('mDT',i))$Ncount, na.rm=TRUE), digits=3))

# Get number of reads with x TC conv (use 30 as arbitrary max)
for (j in c(0:30)){
assign(paste0(paste0('TC',j,'_'),i),nrow(get(paste0('mDT',i))[get(paste0('mDT',i))$TC_mismatches == j]))
}
for (j in c(0:30)){
assign(paste0(paste0(second,j,'_'),i),nrow(get(paste0('mDT',i))[get(paste0('mDT',i))[[mmcolumn]] == j]))
}

# Get total number of reads 
assign(paste0('tot',i),nrow(get(paste0('mDT',i))))

# Get fragment lengths and T counts and G counts
assign(paste0('lengths', i),nchar(get(paste0('mDT',i))$seq))
assign(paste0('TcountsDist', i),get(paste0('mDT',i))$Tcount)
assign(paste0('NcountsDist', i),get(paste0('mDT',i))$Ncount)
}

} # End of first big loop, outside of all if statements


if (task != 'readcountsOnly') {
# Get 2d joint frequency distribution for T counts and TC conversions
for (i in c(1:NumSamps)){
DT <- get(paste0('mDT',i))

JointFreqTC <- table(DT$Tcount,DT$TC_mismatches)
write.table(JointFreqTC, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_TcountANDTCconv_',samples[i],'_otherMMfilt_',stringentFilter,'.txt'), sep=("\t"), quote=FALSE, col.names=NA)

if (second == 'GA') { 
fname = '_GcountANDGAconv_'
} else if (second == 'CT') {
fname = '_CcountANDCTconv_'
}

JointFreqNN <- table(DT$Ncount,DT[[mmcolumn]])
write.table(JointFreqNN, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region, fname, samples[i],'_otherMMfilt_',stringentFilter,'.txt'), sep=("\t"), quote=FALSE, col.names=NA)


# Get list of read names in this subset
# readNames=sapply(strsplit(DT$readName, '/'),"[[",1)
# write.table(readNames, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_readNames_',samples[i],'maxTcount',maxTcount,'.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)

}

if (second == 'GA') { 
f2name = '_LengthAndTandGcountDist'
title = ' G counts distribution'
} else if (second == 'CT') {
f2name = '_LengthAndTandCcountDist'
title = ' C counts distribution'
}

# Make plots for fragment length and T counts
samples = samples # shortnames
pdf(paste0(path, 'MMfrequency/',Exp,'_', MapMethod, '_frag_',region, f2name, '_otherMMfilt_',stringentFilter,'.pdf'), width = 4, height = 12)
par(mfrow=c(3,1), cex.lab = 1)

for (i in c(1:NumSamps)){
hist(get(paste0('lengths',i)), main= paste0(samples[i], ' fragments length distribution'), xlab = paste0(samples[i], ' lengths'), breaks = 50)
assign(paste0('Thist',i), hist(get(paste0('TcountsDist',i)), main= paste0(samples[i], ' T counts distribution'), xlab = paste0('Number of Ts'), breaks = 135))
assign(paste0('Nhist',i), hist(get(paste0('NcountsDist',i)), main= paste0(samples[i], title), xlab = paste0('Number of ',ntF,'s'), breaks = 135))
}
# measure skewness and kurtosis
# library(moments)
# assign(paste0('skew',i), skewness(get(paste0('mDT',i))$Tcount))
# assign(paste0('kurt',i), kurtosis(get(paste0('mDT',i))$Tcount))
dev.off()





# Plot total number of mismatches for TC > n
# n=2
# for (i in c(1:NumSamps)){
# assign(paste0('tot_MM', i), c()) }
# 
# for (i in c(1:NumSamps)){
# assign(paste0('tot_MM', i), (get(paste0('mDT', i))[TC_mismatches>n]$tot_mismatches) - (get(paste0('mDT', i))[TC_mismatches>n]$TC_mismatches))
# }
# 
# pdf(paste0(path, 'MMfrequency/', MapMethod, '_frag_',reads,'_totMMinReadsWithGrThan', n,'TC.pdf'), width = 5, height = 6)
# 
# cols = c('dodgerblue','darkblue', 'violet', 'forestgreen', 'yellowgreen', 'gold1', 'orange', 'red', 'brown', 'grey70', 'grey70', 'grey70')
# cols = cols[1:NumSamps]
# boxplot(tot_MM1,tot_MM2,tot_MM3,tot_MM4,tot_MM5,tot_MM6,tot_MM7,tot_MM8,tot_MM9,tot_MM10,tot_MM11,tot_MM12, ylab = 'Total non-TC mismatches',col = cols, main = paste0('Total mismatches for reads with > ',n,'TC'),xaxt = "n", pch=16, ylim=c(0,15), cex = .5)
# axis(side = 1, at = c(seq(1,NumSamps)), labels = samples, tick = FALSE, las=2)
# 
# dev.off()


# Make fasta to look at reads with >n TC
# library(seqinr)
# n=1
# for (i in c(1:NumSamps)){
# assign(paste0('highTCseq',i),get(paste0('mDT',i))[get(paste0('mDT',i))$TC_mismatches > n]$seq)
# assign(paste0('fastanames',i), seq(1, length(get(paste0('highTCseq',i)))))
# write.fasta(as.list(get(paste0('highTCseq',i))), get(paste0('fastanames',i)), file.out=paste0(path, 'MMfrequency/highTCseqs', samples[i],'.fasta'), open='w')
# }

# Make list of readnames to with >n TC to filter bam (use 'ExtractSelectReadsByNameForGRANDSLAM.sh')

# Get list of read names in this subset
# n=1
# for (i in c(1:NumSamps)) {
# DT <- get(paste0('mDT',i))
# highTCseq <- DT[DT$TC_mismatches > n]$readName
# readNames=sapply(strsplit(highTCseq, '/'),"[[",1)
# write.table(readNames, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_readNames_',samples[i],'GrTh',n,'TCperRead.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
# }
# n=2
# for (i in c(1:NumSamps)) {
# DT <- get(paste0('mDT',i))
# highTCseq <- DT[DT$TC_mismatches > n]$readName
# readNames=sapply(strsplit(highTCseq, '/'),"[[",1)
# write.table(readNames, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_readNames_',samples[i],'GrTh',n,'TCperRead.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
# }
# n=3
# for (i in c(1:NumSamps)) {
# DT <- get(paste0('mDT',i))
# highTCseq <- DT[DT$TC_mismatches > n]$readName
# readNames=sapply(strsplit(highTCseq, '/'),"[[",1)
# write.table(readNames, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_readNames_',samples[i],'GrTh',n,'TCperRead.txt'), sep=("\t"), quote=FALSE, col.names=FALSE, row.names=FALSE)
# }



# Plot average # Ts and Gs per read with std dev and num with TC and GA MM

# Get max number of TC conversions in data
TC_mm = c()
NN_mm = c()
for (i in c(1:NumSamps)) {
TC_mm = c(TC_mm, get(paste0('mDT',i))$TC_mismatches)
NN_mm = c(NN_mm, get(paste0('mDT',i))[[mmcolumn]])
}
maxTC = max(TC_mm)
maxNN = max(NN_mm)


if (second == 'GA') { 
f3name = '_frag_AveTsGs_'
} else if (second == 'CT') {
f3name = '_frag_AveTsCs_'
}

pdf(paste0(path, 'MMfrequency/',Exp,'_', MapMethod, f3name,region,'_otherMMfilt_',stringentFilter,'_NumReadsWithMM.pdf'), width = 20, height = 13)

par(mfrow=c(5,1), cex.lab = 1)

# For barplot T counts rates

Tcountslist=c()
Tsdevlist=c()
Ncountslist=c()
Nsdevlist=c()
totslist=c()

for (i in c(1:NumSamps)){
Tcountslist=c(Tcountslist, get(paste0('meanTcount',i)))
Tsdevlist=c(Tsdevlist, get(paste0('Tstddev',i)))
Ncountslist=c(Ncountslist, get(paste0('meanNcount',i)))
Nsdevlist=c(Nsdevlist, get(paste0('Nstddev',i)))
totslist=c(totslist, get(paste0('tot',i)))
}

Tcounts <- as.numeric(Tcountslist) 
Tsdevs <- as.numeric(Tsdevlist) 
Ncounts <- as.numeric(Ncountslist) 
Nsdevs <- as.numeric(Nsdevlist) 
tots <-  as.numeric(totslist)
minReadCounts = min(totslist[1:(length(totslist))])

####################################################
#### For modifying the distribution of T counts ###
####################################################

# Find sample and with lowest mean T counts and get distribution
# lowestT = min(Tcounts)
# lowestTsample = samples[which(Tcounts == lowestT)]
# index_lowestTsample = which(Tcounts == lowestT)
# lowestTdist <- data.table(table(get(paste0('mDT',index_lowestTsample))$Tcount))
# # Will want to start filling from the highest counts in distribution
# sortdist <- lowestTdist[rev(order(N))]
# # Express all counts in relation to highest
# sortdist$modefrac <- sortdist$N/sortdist[1,]$N
# 
# # Start file to plot new histograms to check new distributions
# pdf(paste0(path, 'MMfrequency/',Exp,'_', MapMethod, '_frag_',region,'_modTcountDist.pdf'), width = 4, height = 4)
# 
# 
# # Make new 2d joint frequency matrix using this distribution
# # Get 2d joint frequency distribution for T counts and TC conversions
# for (i in c(1:NumSamps)){
# DT <- get(paste0('mDT',i))
# # Sum of all reads that will be in distribution
# tot <- nrow(DT[Tcount<=max(lowestTdist$V1)])
# 
# # Make the data table with new distribution
# newDT=data.table()
# modeTnum = as.numeric(sortdist[1,]$V1)
# newDT <- rbind(newDT, DT[Tcount == modeTnum])
# modecount <- nrow(newDT)
# 
# for (k in c(2:nrow(sortdist))) {
# Tnum = as.numeric(sortdist[k,]$V1)
# frac = as.numeric(sortdist[k,]$modefrac)
# newDT <- rbind(newDT, DT[Tcount == Tnum][sample(.N,min(c(nrow(DT[Tcount == Tnum]), round(frac*modecount,0))))])
# }
# 
# 
# # TC
# JointFreqTC <- table(newDT$Tcount,newDT$TC_mismatches)
# write.table(JointFreqTC, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_TcountANDTCconv_modTcountDist_' ,samples[i],'.txt'), sep=("\t"), quote=FALSE, col.names=NA)
# 
# Thist <- hist(newDT$Tcount, main= paste0(samples[i], ' modified T counts distribution'), xlab = paste0('Number of Ts'), breaks = 135)
# 
# }
# dev.off()

##########################



ylimitsCounts = c(0,40)

TNcounts = matrix(c(Tcounts,Ncounts), nrow = 2, byrow=TRUE)
TNsdevs = matrix(c(Tsdevs, Nsdevs), nrow = 2, byrow=TRUE)

xx = barplot(TNcounts, beside=TRUE, names.arg=samples, col=c('red', 'orange'), cex.names = .8, ylab = paste0('Mean number T/',ntF,' in fragment (std dev)'), main = paste0('T/',ntF,' counts per fragment'), ylim=ylimitsCounts) # 

# Add std dev as text to top of bars
text(x = xx, y = TNcounts, label = TNsdevs, pos = 3, cex = 0.8, col = "black")
# Add total number of reads as text to bottom of bars
text(x = xx[1,]+.5, y = 1, label = tots, pos = 3, cex = 0.8, col = "black")
# Add mean value near top of bar
text(x = xx, y = TNcounts-5, label = TNcounts, pos = 3, cex = 0.8, col = "black")


for (j in c(0:maxTC)){
assign(paste0('TClist',j),c())
assign(paste0('NNlist',j),c())
}

for (j in c(0:maxTC)){
for (i in c(1:NumSamps)){
assign(paste0('TClist',j),c(get(paste0('TClist',j)), get(paste0(paste0('TC',j,'_'),i))))
assign(paste0('NNlist',j),c(get(paste0('NNlist',j)), get(paste0(paste0(second,j,'_'),i))))
}
}

# For frequency of each count
for (j in c(0:maxTC)){
assign(paste0('TC',j,'s'),as.numeric(get(paste0('TClist',j))))
assign(paste0(second,j,'s'),as.numeric(get(paste0('NNlist',j))))
}

Tvec = c()
Nvec = c()
for (j in c(0:maxTC)){
Tvec = c(Tvec,get(paste0('TC',j,'s')))
Nvec = c(Nvec,get(paste0(second,j,'s')))
}
TCs <- matrix(Tvec, nrow=NumSamps, ncol=maxTC+1)
NNs <- matrix(Nvec, nrow=NumSamps, ncol=maxTC+1)


if (NumSamps > 7) {
cols = c('dodgerblue','darkblue', 'blue','darkorchid','gold1', 'forestgreen','green', 'yellowgreen', 'yellow','violet', 'orange', 'red', 'red3','brown')} # skyblue
if (NumSamps == 6 | NumSamps == 7) {
cols = c('dodgerblue','violet', 'forestgreen', 'gold1', 'orange', 'red', 'brown')}
if (NumSamps == 9) {
cols = c('dodgerblue','darkblue','violet', 'forestgreen','green', 'yellow', 'orange', 'red', 'brown')} # skyblue
if (NumSamps <8) {
cols = c('dodgerblue','violet', 'forestgreen', 'gold1', 'orange', 'red', 'brown')}
if (NumSamps == 5) {
cols = c('dodgerblue', 'forestgreen', 'gold1', 'orange', 'red')}
if (NumSamps == 4) {
cols = c('dodgerblue', 'green', 'gold1', 'red')}
if (Exp == 'TL5' & NumSamps == 8) {
cols = c('dodgerblue', 'green', 'gold1', 'red','dodgerblue', 'green', 'gold1', 'red')} # skyblue
if (Exp == 'TL6') {
cols = c('grey60', 'yellow', 'dodgerblue', 'forestgreen')}
if (Exp == 'TL13') {
cols = rep(c('dodgerblue','violet', 'forestgreen', 'gold1', 'orange', 'red'), 4)}

cols = cols[1:NumSamps]


ylimitsBar1 = c(0, 1.2*max(TC1s))

# Plot for TC (unnormalized)
names = as.character(c(0:maxTC))
xxx = barplot(TCs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, ylab = 'Frequency', main = 'Number of fragments with n TC conversions', legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBar1) 

text(x = xxx+.3, y = TCs+(max(ylimitsBar1/50)), label = TCs, pos = 3, cex = 0.6, srt=90,col = 'black')

# Plot for NN (unnormalized)
xxx = barplot(NNs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, ylab = 'Frequency', main = paste0('Number of fragments with n ',second,' conversions'), legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBar1) 

text(x = xxx+.3, y = NNs+(max(ylimitsBar1/50)), label = NNs, pos = 3, cex = 0.6, srt=90,col = 'black')



# For frequency of each count, normalized to per min reads
# Make matrix with readcounts to divide TCs matrix
readcountMat = matrix(tots, nrow=NumSamps, ncol=maxTC+1)
NormTCs = round(TCs/readcountMat*minReadCounts, 0)
NormNNs = round(NNs/readcountMat*minReadCounts, 0)

ylimitsBarNorm = c(0, 1.1*max(NormTCs[,2]))

# Plot for TC (normalized)
xxx = barplot(NormTCs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, , ylab = 'Frequency', main = paste0('Number of fragments with n TC conversions\n normalized per ',minReadCounts,' reads'), legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBarNorm) 

text(x = xxx+.3, y = NormTCs+(max(ylimitsBarNorm)/50), label = NormTCs, pos = 3, cex = 0.6, srt=90,col = 'black')

# Plot for NN (normalized)
xxx = barplot(NormNNs, beside=TRUE, names.arg=names, col=cols, cex.names = 1, , ylab = 'Frequency', main = paste0('Number of fragments with n ',second,' conversions\n normalized per ',minReadCounts,' reads'), legend.text=samples, args.legend = list(x='topright', bty = 'n', fill = cols, cex = .8),ylim=ylimitsBarNorm) 

text(x = xxx+.3, y = NormNNs+(max(ylimitsBarNorm)/50), label = NormNNs, pos = 3, cex = 0.6, srt=90,col = 'black')

dev.off()

# Write conversion frequencies to table
TCconvFreq <- data.table(t(TCs))
colnames(TCconvFreq) <- samples
TCconvFreq[, TCperRead := seq(0,maxTC)]
setcolorder(TCconvFreq, c('TCperRead',samples))
# colnames(TCconvFreq) <- samples
write.table(TCconvFreq, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region,'_otherMMfilt_',stringentFilter,'_TCconvFreq.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)

NNconvFreq <- data.table(t(NNs))
colnames(NNconvFreq) <- samples
NNconvFreq[, NNperRead := seq(0,maxTC)]
setcolorder(NNconvFreq, c('NNperRead',samples))
# colnames(TCconvFreq) <- samples

if (second == 'GA') { 
f4name = '_GAconvFreq'
} else if (second == 'CT') {
f4name = '_CTconvFreq'
}
write.table(NNconvFreq, file=paste0(path, 'MMfrequency/',Exp,'_',MapMethod,'_frag_',region, f4name,'_otherMMfilt_',stringentFilter,'.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)
}
# source('/n/groups/churchman/mc348/TimelapseSeq/Scripts/MismatchFrequencyTCandGA.R')
