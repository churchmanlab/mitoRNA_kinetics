# First run the following on .gtf to choose genes and get just the name and coordinate columns for transcripts:
# infile="ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes.gtf"
# outfile="MitoGenes_transcripts.txt"
# grep -wFf /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Mito\ genes_list_newSymbols.txt $infile | awk -F'[\t|"]' '{print $1,$2,$3,$4,$5,$18,$22}' | grep 'transcript' |grep 'protein_coding' | awk '{print $6,$4,$5}' > $outfile


library(data.table)

path <- '/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/'

infile <- 'MitoGenes_transcripts.txt' 

outfile <- 'MitoGenes_transcripts_lengths.txt' 

DT <- data.table(read.table(paste0(path, infile), sep = ' ', header=FALSE))

DT <- DT[ , Length := V3-V2]

DT <- DT[!grep('-',DT$V1)]


> summary(DT$Length)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    249    6000   12014   20180   22438  132876 
