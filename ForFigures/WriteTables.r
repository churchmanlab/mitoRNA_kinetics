# Abundance

type = 'RPKM' # reads, RPK, RPKM, RPKS 

TL3 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL3_2020_09/RPK/TL3_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

TL4 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL4_combined/RPK/TL4_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

TL1 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Hela_TL1_2020_02/RPK/TL1_t5MMinformed4_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

TL10 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL10_LRPPRC_2022_07/RPK/TL10_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

TL12 <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/HEK_TL12_LRPPRC_2022_08/RPK/TL12_t5MTMMinformed6_featureCounts_multi_',type,'_noPseudo.txt'), sep="\t",header=TRUE, stringsAsFactors=FALSE))

DTs = c('TL3','TL4','TL1','TL10','TL10','TL12', 'TL12')
samp = c('WT','WT','WT','WT','LRP','WT','LRP')

mitogenes = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6')

mitolist = c('MT-ND1',  'MT-ND2',  'MT-ND3', 'MT-ND4L-4', 'MT-ND5', 'MT-CYB', 'MT-CO1',  'MT-CO2',  'MT-CO3', 'MT-ATP8-6','MT-ND6', 'MT-RNR2', TL3[grep('MT-anti', DT$GeneName)]$GeneName, '7S') # this one includes RNR2 and anti genes

nucgenes = data.table(read.table('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/Mito genes_newSymbols.txt', sep="\t",header=TRUE, stringsAsFactors=FALSE))$Mito.genes

for (i in c(1:length(DTs))) {

DT = get(DTs[[i]])

column=2

if (samp[i] == 'WT') {
Mito = DT[GeneName %in% mitolist][,c(1,2)]
Nuc = DT[GeneName %in% nucgenes][,c(1,2)]
}
if (samp[i] == 'LRP') {
Mito = DT[GeneName %in% mitolist][,c(1,3)]
Nuc = DT[GeneName %in% nucgenes][,c(1,3)]
}
# Remove 0 values
# Nuc_1 <- Nuc_1[Nuc_1[[1]]!=0]



# For writing files

write.table(setorder(Mito, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/',DTs[i],'_',samp[i],'_mito_Ab.txt'), sep='\t',row.names=F, quote=F)
write.table(setorder(Nuc, cols='GeneName'), paste0('/Users/Mary/Desktop/Data/TimelapseSeq/Comparisons/Values/',DTs[i],'_',samp[i],'_nuc_Ab.txt'), sep='\t', row.names=F, quote=F)

}


# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ForFigures/WriteTables.r')

