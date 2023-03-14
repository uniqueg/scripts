### A. PRE-REQUISITES
# Initialize and get command line arguments
args <- commandArgs(trailingOnly=TRUE)
## Pass arguments
inp <- "prob_cut_ind_vs_unind_alpha_1"
out <- "prob_cut_ind_vs_unind_alpha_1_log.pdf"
###

### B. LOAD/PREPARE DATA
# Read input matrix
df <- read.table(inp, header=FALSE, row.names=NULL, sep=" ")
# Remove last column resulting from trailing <space> introduced by C++ program "psamepdiffmodel" (by Mihaela Zavolan)
df <- df[,-101]
###

### C. PLOT
# Initiate pdf graphic device
pdf(out, height=6, width=30)
# Generate boxplot
boxplot(df, log="y", names=1:100, cex.axis=0.75, outline=FALSE)
# Close pdf graphic device
dev.off()
###
