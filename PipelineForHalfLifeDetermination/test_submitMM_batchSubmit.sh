# 2019-09-26

# To run after submit_parallel_findMM.sh to catch any temp files that did not get MM files

# use: /n/groups/churchman/mc348/TimelapseSeq/Scripts/test_submitMM.sh 

prefix=$1

# set path
# path=`pwd`/${prefix}/

# remove the scripts for each tmp file
rm -f ${prefix}/*_mm.sh

# list all tmp files to run through
tmp_files=$(GLOBIGNORE="${prefix}/*.bed"; ls ${prefix}/tmp*)

unset GLOBIGNORE

# if _MM.sam file does not exist for the tmp file, print the tmp file
for x in $tmp_files 

do

[ ! -f "${x}_MM.bed" ] && echo "$x"

done


