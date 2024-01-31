# Install

source ~/Python_Venv/bin/activate
unset PYTHONPATH # to keep it from using the old version of numpy
# pip install pyGenomeTracks


# Downsample bam
# sbatch -p short -t 0-02:00 -c 4 --wrap="samtools view -@ 3 -bs 42.05 H_1_6_chrMall_noDups.bam > H_1_6_chrMall_noDups_ds.bam"

# Remove rRNA
sbatch -p short -t 0-02:00 -c 4 --wrap="samtools view -@ 3 -b H_1_6_chrMall_noDups.bam chrM:3305 -o H_1_6_chrMall_noDups_norRNA.bam"


# Convert bam to bed
sbatch -p short -t 0-02:00 --wrap="bamToBed -cigar -i H_1_6_chrMall_noDups_norRNA.bam > H_1_6_chrMall_noDups_norRNA.bed"

# Remove split reads $7 is cigar string
awk '{if ($7 !~ /N/) print $0;}' H_1_6_chrMall_noDups_norRNA.bed > H_1_6_chrMall_noDups_norRNA_noN.bed
awk '{if ($6 == "+") print $0;}' H_1_6_chrMall_noDups_norRNA_noN.bed > H_1_6_chrMall_noDups_norRNA_noN_P.bed


# Just run this once for each sample to make file
# sbatch -p short -t 0-02:00 --wrap="make_tracks_file --trackFiles H_1_6_chrMall_noDups_norRNA_noN_P_sub_ds.bed /n/groups/churchman/hMitoRP/SequenceFiles/chrM_txptends.gtf -o tracks.ini"

# COX1 5'
let s1=5901-30
let s2=5901+3
let s=5880
let e=5930

# COX3 5'
let s1=9206-30
let s2=9206+3
let s=9206-20
let e=9206+20

# ND1 5'
let s1=3307-24 # (there is a tRNA read up to this point)
let s2=3307+3
let s=3307-20
let e=3307+20

# CYB 3'
let s1=15887-35
let s2=15887-3
let s=15887-20
let e=15887+20

# CO2 3'
let s1=8294-35
let s2=8294-3
let s=8294-20
let e=8294+20

# Filter to just reads in range of interest
awk -v s1=$s1 -v s2=$s2 '{if ($2 > s1 && $2 < s2) print $0;}' H_1_6_chrMall_noDups_norRNA_noN_P.bed > H_1_6_chrMall_noDups_norRNA_noN_P_sub.bed
# Downsample to N lines
N=100
shuf -n $N H_1_6_chrMall_noDups_norRNA_noN_P_sub.bed > H_1_6_chrMall_noDups_norRNA_noN_P_sub_ds.bed


# CO1
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:${s}-${e} -o H_1_6_CO1_5p.pdf"

# CO3
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:${s}-${e} -o H_1_6_CO3_5p.pdf"

# ND1
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:${s}-${e} -o H_1_6_ND1_5p.pdf"


# CYB 3'
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:${s}-${e} -o H_1_6_CYB_3p.pdf"

# CO2 3'
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:${s}-${e} -o H_1_6_CO2_3p.pdf"



#ND1
sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:5880-5940 --width 30 -o H_1_6_CO1_5p.pdf"

sbatch -p short -t 0-04:00 --mem=15G --wrap="pyGenomeTracks --tracks tracks.ini --region chrM:5870-5940 -o H_1_6_CO1_5p.pdf"