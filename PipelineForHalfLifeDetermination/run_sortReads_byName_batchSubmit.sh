#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

module load gcc/6.2.0
module load R/4.0.1

prefix=$1
echo $prefix

path=`pwd`/${prefix}

/n/groups/churchman/mc348/TimelapseSeq/Scripts/sortReads_byName.R reads.MM.bed $path


