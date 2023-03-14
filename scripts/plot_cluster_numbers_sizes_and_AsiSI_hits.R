# Load data
#ind_no_size <- read.table("/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/clusters/output/cluster_counts_induced", sep="\t", fill=TRUE)
#unind_no_size <- read.table("/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/clusters/output/cluster_counts_uninduced", sep="\t", fill=TRUE)
#ind_hits <- read.table("/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/output/cluster_AsiSI_overlaps_induced", sep="\t", fill=TRUE)
#unind_hits <- read.table("/import/bc2/home/zavolan/kanitz/DDRNA/chip-seq/potential_cleavage_sites/output/cluster_AsiSI_overlaps_uninduced", sep="\t", fill=TRUE)
ind_no_size <- read.table("cluster_counts_induced", sep="\t", fill=TRUE)
unind_no_size <- read.table("cluster_counts_uninduced", sep="\t", fill=TRUE)
ind_hits <- read.table("cluster_AsiSI_overlaps_induced", sep="\t", fill=TRUE)
unind_hits <- read.table("cluster_AsiSI_overlaps_uninduced", sep="\t", fill=TRUE)

# Sort by window size
ind_no_size_sort <- ind_no_size[order(ind_no_size[,1]), ]
unind_no_size_sort <- unind_no_size[order(unind_no_size[,1]), ]
ind_hits <- ind_hits[order(ind_hits[,1]), ]
unind_hits <- unind_hits[order(unind_hits[,1]), ]

# strange small numbers inside!
ind_hits <- ind_hits[308:3307,]
unind_hits <- unind_hits[12:3011,]

# Subset data to plot
x <- ind_no_size_sort[ ,1]
no_ind <- ind_no_size_sort[ ,2]
no_unind <- unind_no_size_sort[ ,2]
avg_size_ind <- as.integer(rowMeans(ind_no_size_sort[,3:length(ind_no_size_sort)], na.rm=TRUE))
avg_size_unind <- as.integer(rowMeans(unind_no_size_sort[,3:length(unind_no_size_sort)], na.rm=TRUE))
med_size_ind <- as.integer(apply(ind_no_size_sort[,3:length(ind_no_size_sort)], 1, function(x) median(x, na.rm=TRUE)))
med_size_unind <- as.integer(apply(unind_no_size_sort[,3:length(unind_no_size_sort)], 1, function(x) median(x, na.rm=TRUE)))
hits_all_ind <- ind_hits[ ,2]
hits_all_unind <- unind_hits[ ,2]
hits_1_ind <- as.integer(apply(ind_hits[,3:length(ind_hits)], 1, function(x) { 
	tab <- table(x) 
	return(as.integer(tab[names(tab) == 1]))
}))
hits_1_unind <- as.integer(apply(unind_hits[,3:length(unind_hits)], 1, function(x) {
	tab <- table(x) 
	return(as.integer(tab[names(tab) == 1]))
}))
size_dist_10000_ind <- na.omit(as.integer(ind_no_size_sort[10,3:length(ind_no_size_sort[10, ])]))
size_dist_10000_unind <- na.omit(as.integer(unind_no_size_sort[10,3:length(unind_no_size_sort[10, ])]))
size_dist_100000_ind <- na.omit(as.integer(ind_no_size_sort[100,3:length(ind_no_size_sort[100, ])]))
size_dist_100000_unind <- na.omit(as.integer(unind_no_size_sort[100,3:length(unind_no_size_sort[100, ])]))
size_dist_1000000_ind <- na.omit(as.integer(ind_no_size_sort[1000,3:length(ind_no_size_sort[1000, ])]))
size_dist_1000000_unind <- na.omit(as.integer(unind_no_size_sort[1000,3:length(unind_no_size_sort[1000, ])]))

# Log scale
x_log <- log(x,10)
x_log2 <- log(x,2)
no_ind_log <- log(no_ind,10)
no_unind_log <- log(no_unind,10)
avg_size_ind_log <- log(avg_size_ind,10)
avg_size_unind_log <- log(avg_size_unind,10)
med_size_ind_log <- log(med_size_ind,10)
med_size_unind_log <- log(med_size_unind,10)
hits_all_ind_log2 <- log(hits_all_ind,2)
hits_all_unind_log2 <- log(hits_all_unind,2)
hits_1_ind_log2 <- log(hits_1_ind,2)
hits_1_unind_log2 <- log(hits_1_unind,2)
size_dist_10000_ind_log <- log(size_dist_10000_ind, 10)
size_dist_10000_unind_log <- log(size_dist_10000_unind, 10)
size_dist_100000_ind_log <- log(size_dist_100000_ind, 10)
size_dist_100000_unind_log <- log(size_dist_100000_unind, 10)
size_dist_1000000_ind_log <- log(size_dist_1000000_ind, 10)
size_dist_1000000_unind_log <- log(size_dist_1000000_unind, 10)



# Plot number of clusters per window size
pdf("cluster_counts_per_window_size.pdf", width = 8, height = 8)
plot(x_log, no_ind, type="l", col="blue", main="Number of clusters per window size, induced vs uninduced", xlab="Window size", ylab="Number of clusters")
lines(x_log, no_unind, col="brown")
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

# Plot number of clusters with AsiSI sites
pdf("clusters_with_AsiSI_sites_per_window_size.pdf", width = 8, height = 8)
plot(x_log, hits_all_ind, type="l", col="blue", main="Number of clusters with AsiSI sites per window size, induced vs uninduced", xlab="Window size", ylab="Number of clusters with AsiSI sites")
lines(x_log, hits_all_unind, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot proportions of clusters with AsiSI sites
pdf("clusters_with_AsiSI_sites_per_window_size_proportions.pdf", width = 8, height = 8)
plot(x_log, hits_all_ind/no_ind, type="l", col="blue", main="Proportions of clusters with AsiSI sites per window size, induced vs uninduced", xlab="Window size", ylab="Proportion of clusters with AsiSI sites")
lines(x_log, hits_all_unind/no_unind, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot number of clusters with exactly one AsiSI sites
pdf("clusters_with_1_AsiSI_site_per_window_size.pdf", width = 8, height = 8)
plot(x_log, hits_1_ind, type="l", col="blue", main="Number of clusters with one AsiSI site per window size, induced vs uninduced", xlab="Window size", ylab="Number of clusters with one AsiSI site")
lines(x_log, hits_1_unind, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot proportions of clusters with exactly one AsiSI sites
pdf("clusters_with_1_AsiSI_site_per_window_size_proportions.pdf", width = 8, height = 8)
plot(x_log, hits_1_ind/no_ind, type="l", col="blue", main="Proportions of clusters with one AsiSI site per window size, induced vs uninduced", xlab="Window size", ylab="Proportion of clusters with one AsiSI site")
lines(x_log, hits_1_unind/no_unind, col="brown")
legend("topleft", legend=c("induced", "uninduced"), lwd=1, col=c("blue", "brown"))
dev.off()

# Plot number of clusters with any and exactly one AsiSI site
pdf("clusters_with_any_or_1_AsiSI_site_per_window_size.pdf", width = 8, height = 8)
plot(x_log, hits_all_ind, type="l", col="blue", main="Number of clusters with any or one AsiSI site per window size, induced vs uninduced", xlab="Window size", ylab="Number of clusters with AsiSI sites")
lines(x_log, hits_all_unind, col="brown")
lines(x_log, hits_1_ind, col="lightblue")
lines(x_log, hits_1_unind, col="brown3")
legend("topleft", legend=c("Any number of AsiSI sites, induced", "Any number of AsiSI sites, uninduced", "1 AsiSI site, induced", "1 AsiSI site, uninduced"), lwd=1, col=c("blue", "brown3", "lightblue", "orange"))
dev.off()

# Plot proportions of clusters with any and exactly one AsiSI site
pdf("clusters_with_any_or_1_AsiSI_site_per_window_size_proportions.pdf", width = 8, height = 8)
plot(x_log, hits_all_ind/no_ind, type="l", col="blue", main="Proportions of clusters with any or one AsiSI site per window size, induced vs uninduced", xlab="Window size", ylab="Proportion of clusters with AsiSI sites")
lines(x_log, hits_all_unind/no_unind, col="brown")
lines(x_log, hits_1_ind/no_ind, col="lightblue")
lines(x_log, hits_1_unind/no_unind, col="brown3")
legend("topleft", legend=c("Any number of AsiSI sites, induced", "Any number of AsiSI sites, uninduced", "1 AsiSI site, induced", "1 AsiSI site, uninduced"), lwd=1, col=c("blue", "brown3", "lightblue", "orange"))
dev.off()

# Plot size distribution for 10^4
pdf("cluster_size_distribution_window_size_10000_induced.pdf", width = 12, height = 6)
hist(size_dist_10000_ind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^4, induced", xlab="Cluster size", ylim=c(0,1))
dev.off()
pdf("cluster_size_distribution_window_size_10000_uninduced.pdf", width = 12, height = 6)
hist(size_dist_10000_unind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^4, uninduced", xlab="Cluster size", ylim=c(0,1))
dev.off()

# Plot size distribution for 10^5
pdf("cluster_size_distribution_window_size_100000_induced.pdf", width = 12, height = 6)
hist(size_dist_100000_ind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^5, induced", xlab="Cluster size", ylim=c(0,1))
dev.off()
pdf("cluster_size_distribution_window_size_100000_uninduced.pdf", width = 12, height = 6)
hist(size_dist_100000_unind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^5, uninduced", xlab="Cluster size", xlim=c(3.5,8), ylim=c(0,1))
dev.off()

# Plot size distribution for 10^6
pdf("cluster_size_distribution_window_size_1000000_induced.pdf", width = 12, height = 6)
hist(size_dist_1000000_ind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^6, induced", xlab="Cluster size", ylim=c(0,1))
dev.off()
pdf("cluster_size_distribution_window_size_1000000_uninduced.pdf", width = 12, height = 6)
hist(size_dist_1000000_unind_log, probability=TRUE, breaks=250, main="Cluster size distribution for window size 10^6, uninduced", xlab="Cluster size", ylim=c(0,1))
dev.off()