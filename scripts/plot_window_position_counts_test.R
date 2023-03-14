ind_yH2Ax_ind <- read.table("AsiSI_in_clusters_100000_induced_window_positon_counts_yH2Ax_induced", row.names=1, fill=TRUE)
ind_chr_ind <- read.table("AsiSI_in_clusters_100000_induced_window_positon_counts_chromatin_induced", row.names=1, fill=TRUE)
unind_yH2Ax_unind <- read.table("AsiSI_in_clusters_100000_uninduced_window_positon_counts_yH2Ax_uninduced", row.names=1, fill=TRUE)
unind_chr_unind <- read.table("AsiSI_in_clusters_100000_uninduced_window_positon_counts_chromatin_uninduced", row.names=1, fill=TRUE)

win_size_l <- 100000
win_size_s <- 10000


df <- ind_yH2Ax_ind

### Plots for indicated window size
## Subset data
# Generate logical vector indicating which columns to subset
cols_l <- c	(rep( FALSE, ( (ncol(df) - win_size_l) / 2 )), rep( TRUE , (win_size_l + 1)), rep( FALSE, (((ncol(df) - win_size_l) / 2))))
# Subset dataframe
ind_yH2Ax_ind_l <- subset(ind_yH2Ax_ind, select=cols_l)
ind_chr_ind_l <- subset(ind_chr_ind, select=cols_l)
unind_yH2Ax_unind_l <- subset(unind_yH2Ax_unind, select=cols_l)
unind_chr_unind_l <- subset(unind_chr_unind, select=cols_l)
# Generate x axis labels (window positions centered around 0)
x_l <- seq(-(win_size_l / 2), win_size_l / 2, 1)
##

## Plot average window counts
# Compute averages
ind_yH2Ax_ind_l <- colSums(ind_yH2Ax_ind_l) / nrow(ind_yH2Ax_ind_l)
ind_chr_ind_l <- colSums(ind_chr_ind_l) / nrow(ind_chr_ind_l)
unind_yH2Ax_unind_l <- colSums(unind_yH2Ax_unind_l) / nrow(unind_yH2Ax_unind_l)
unind_chr_unind_l <- colSums(unind_chr_unind_l) / nrow(unind_chr_unind_l)

# Compute differences
induced <- ind_yH2Ax_ind_l - ind_chr_ind_l
uninduced <- unind_yH2Ax_unind_l - unind_chr_unind_l

# Open pdf printer device
pdf("environment_100000_induced.pdf", width = 24, height = 8)
# Plot
plot(x_l, induced, type="l", main="Induced", xlab="Position relative to AsiSI motifs", ylab="Read count average (FG - BG reads)")
# Close pdf printer
dev.off()
##

# Open pdf printer device
pdf("environment_100000_uninduced.pdf", width = 24, height = 8)
# Plot
plot(x_l, uninduced, type="l", main="Uninduced", xlab="Position relative to AsiSI motifs", ylab="Read count average (FG - BG reads)")
# Close pdf printer
dev.off()
##


### Plots for innermost 10% of indicated window size
## Subset data
# Generate logical vector indicating which columns to subset
cols_s <- c(rep(FALSE, ((win_size_l - win_size_s) / 2)), rep(TRUE, (win_size_s + 1)), rep(FALSE, ((win_size_l - win_size_s) / 2)))
# Subset dataframe
induced_s <- subset(induced, subset=cols_s)
uninduced_s <- subset(uninduced, subset=cols_s)
# Generate x axis labels (window positions centered around 0)
x_s <- seq(-(win_size_s / 2), (win_size_s / 2), 1) 
##

# Open pdf printer device
pdf("environment_10000_induced.pdf", width = 24, height = 8)
# Plot
plot(x_s, induced_s, type="l", main="Induced", xlab="Position relative to AsiSI motifs", ylab="Read count average (FG - BG reads)")
# Close pdf printer
dev.off()
##

# Open pdf printer device
pdf("environment_10000_uninduced.pdf", width = 24, height = 8)
# Plot
plot(x_s, uninduced_s, type="l", main="Uninduced", xlab="Position relative to AsiSI motifs", ylab="Read count average (FG - BG reads)")
# Close pdf printer
dev.off()
##
