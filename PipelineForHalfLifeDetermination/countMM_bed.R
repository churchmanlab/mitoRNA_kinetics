#!/n/app/R/4.0.1/bin/Rscript

# Count occurances of all types of mismatches within SLAM-seq data, using output from findMismatches_complete.R
# For Brendan, written on 3/21/19 by K. Lachance

# Modified by Brendan on 11/01/2019 to work with bed files generated by findMismatches_complete_19_10_23.R

# Use: ./countMM.R *_MM.sam

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
path <- paste0(args[2],"/")

# Read in dataframe
dat <- read.table(paste0(path,dat_file), sep="\t", stringsAsFactors = FALSE)

# Create a dataframe to store output (4x4 matrix with rows for each possible genomic nt and columns for each possible read nt)
nt = c("A", "C", "T", "G")
out <- data.frame(matrix(0, nrow = 4, ncol = 4), stringsAsFactors = FALSE)
row.names(out) <- nt # Genome nt
colnames(out) <- nt # Read nt

# Loop through each line (read)
for (i in 1:nrow(dat)) {

    # Check if there are any mismatches
    if(dat[i,18] > 0) {
    
	# Get read and genome nucleotides
	read.nt = dat[i, 13]
	gen.nt = dat[i, 15]

	# Split by commas
	r.nt = gsub(" ", "", unlist(strsplit(read.nt, ",")), fixed = TRUE)
	g.nt = toupper(gsub(" ", "", unlist(strsplit(gen.nt,	",")), fixed = TRUE))

	# Skip any "N" nucleotides
	N = r.nt != "N"
	r.nt = r.nt[N]
	g.nt = g.nt[N]

	# Add to output dataframe based on mismatch 
	for (n in 1:length(r.nt)) {
	    
	    r = r.nt[n]
	    g = g.nt[n]

#	    if (r == g) {
#	       print(dat[i,]) # Print if the MM is the same (called incorrectly)
#	    }

	    out[g, r] = out[g, r] + 1

	} # End of for loop

    } # End of if statement

} # End of for loop

# Report results
cat(paste("CA=", out["C", "A"], "\n", sep=""))
cat(paste("CG=", out["C", "G"],	"\n", sep=""))
cat(paste("CT=", out["C", "T"],	"\n", sep=""))
cat(paste("CC=", out["C", "C"],	"\n", sep=""))
cat(paste("AC=", out["A", "C"],	"\n", sep=""))
cat(paste("AG=", out["A", "G"],	"\n", sep=""))
cat(paste("AT=", out["A", "T"],	"\n", sep=""))
cat(paste("AA=", out["A", "A"],	"\n", sep=""))
cat(paste("GC=", out["G", "C"],	"\n", sep=""))
cat(paste("GA=", out["G", "A"],	"\n", sep=""))
cat(paste("GT=", out["G", "T"],	"\n", sep=""))
cat(paste("GG=", out["G", "G"],	"\n", sep=""))
cat(paste("TC=", out["T", "C"],	"\n", sep=""))
cat(paste("TA=", out["T", "A"],	"\n", sep=""))
cat(paste("TG=", out["T", "G"],	"\n", sep=""))
cat(paste("TT=", out["T", "T"],	"\n", sep=""))