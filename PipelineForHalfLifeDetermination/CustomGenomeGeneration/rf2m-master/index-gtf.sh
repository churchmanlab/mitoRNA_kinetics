#!/bin/bash                                                                                                                                    

#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###

if [ "$#" -ne 1 ]; then
  echo -e "Incorrect number of parameters! Usage:\n    index-gtf.sh <file.gtf(.gz)>" >&2
  exit 1
fi

gtf="$1"

if [[ $gtf =~ \.gz$ ]]; then
  output=${gtf%.gtf.gz}.sorted.gtf.gz
  zcat $gtf | awk '!( $0 ~ /^#/ )' | sort --parallel=4 -S4G -k1,1 -k4,4n | /n/groups/churchman/mc348/programs/tabix-0.2.6/bgzip -c > $output
else
  output=${gtf%.gtf}.sorted.gtf.gz
  cat $gtf | awk '!( $0 ~ /^#/ )' | sort --parallel=4 -S4G -k1,1 -k4,4n | /n/groups/churchman/mc348/programs/tabix-0.2.6/bgzip -c > $output
fi
/n/groups/churchman/mc348/programs/tabix-0.2.6/tabix $output

