# Load data
ind <- read.table("cluster_counts_induced", sep="\t", fill=TRUE)
unind <- read.table("cluster_counts_uninduced", sep="\t", fill=TRUE)

# Sort by window size
ind_sort <- ind[ order(ind[,1]), ]
unind_sort <- unind[ order(unind[,1]), ]

# Subset data to plot
x <- ind_sort[ ,1]
no_ind <- ind_sort[ ,2]
no_unind <- unind_sort[ ,2]
avg_size_ind <- as.integer(rowMeans(ind_sort[,3:length(ind_sort)], na.rm=TRUE))
avg_size_unind <- as.integer(rowMeans(unind_sort[,3:length(unind_sort)], na.rm=TRUE))
med_size_ind <- as.integer(apply(ind_sort[,3:length(ind_sort)], 1, function(x) median(x, na.rm=TRUE)))
med_size_unind <- as.integer(apply(unind_sort[,3:length(unind_sort)], 1, function(x) median(x, na.rm=TRUE)))

# Log scale
x_log <- log(x,10)
no_ind_log <- log(no_ind,10)
no_unind_log <- log(no_unind,10)
avg_size_ind_log <- log(avg_size_ind,10)
avg_size_unind_log <- log(avg_size_unind,10)
med_size_ind_log <- log(med_size_ind,10)
med_size_unind_log <- log(med_size_unind,10)

# Plot number of clusters per window size
pdf("cluster_counts_per_window_size.pdf", width = 8, height = 8)
plot(x_log, no_ind_log, type="l", col="blue", main="Number of clusters per window size, induced vs uninduced", xlab="Window size", ylab="Number of clusters")
lines(x_log, no_unind_log, col="brown")
legend("topright", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot average size per window size
pdf("cluster_size_average_per_window_size.pdf", width = 8, height = 8)
plot(x_log, avg_size_ind_log, type="l", col="blue", main="Average cluster size per window size, induced vs uninduced", xlab="Window size", ylab="Average cluster size")
lines(x_log, avg_size_unind_log, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot median size per window size
pdf("cluster_size_median_per_window_size.pdf", width = 8, height = 8)
plot(x_log, avg_size_ind_log, type="l", col="blue", main="Median cluster size per window size, induced vs uninduced", xlab="Window size", ylab="Median cluster size")
lines(x_log, avg_size_unind_log, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()
