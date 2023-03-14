### Author: Alexander Kanitz
### Created: 20-MAY-2013
### Modified: 20-MAY-2013
### Description: Sums up counts of a position matrix (e.g. resulting from overlap_start_site_counts_window_around_pos_vs_bed.pl and related) in bins of an indicated size
### Arguments: 1. Position matrix [FILE; space-separated]; 2. Window/bin size [INTEGER]; 6. Output file [FILE]
### Output: Bin matrix 
### Usage: Rscript ./sum_win_pos_counts.R ./in_matrix 100 ./out_matrix 

### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
inp <- as.character(args[1])
win <- as.numeric(args[2])
out <- as.character(args[3])
###

### B. READ DATA
# Read input matrix
df <- read.table(inp, row.names=1, fill=TRUE, header=FALSE)
###

### C. CALCULATE BIN POSITIONS
# Extract matrix width
width <- ncol(df)
# Calculate the number of bases that cannot fit in the region with the specified window size
rest <- width %% win
## Calculate offsets for window positions on either side of the region to account for bases that do not fit in windows
offset_start <- floor(rest / 2)
offset_end <- rest - offset_start
## Calculate start and end positions of bins
start <- seq(1 + offset_start, width - offset_end, win)
end <- seq(win + offset_start, width - offset_end, win)
###

### D. SUM COUNTS
# Initiate list to keep the bin counts
red_df <- list()
## Iterate over all bin positions...
for (i in 1:length(start)) {
	# ...and calculate counts
	red_df[[i]] <- rowSums(df[,start[i]:end[i]])
}
# Generate data frame from list of counts
red_df <- do.call(cbind, red_df)
###

### E. WRITE OUTPUT
# Write output matrix to file
write.table(red_df, file=out, quote=FALSE, col.names=FALSE)
###