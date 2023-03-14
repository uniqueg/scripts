### DATA PREPARATION

## INDUCED
# Load data 
i <- read.table("chip_peaks_induced_cont", col.names = "position")
# Subset data
i <- i[,1]
# Plot
pdf("i.pdf", width = 24, height = 8)
hist(i, breaks=30000, main="ChIP peaks across genome, induced", xlab="Continuous Genome Position")
dev.off()

## UNINDUCED
# Load data
u <- read.table("chip_peaks_uninduced_cont", col.names = "position")
# Subset data
u <- u[,1]
# Plot
pdf("u.pdf", width = 24, height = 8)
hist(u, breaks=30000, main="ChIP peaks across genome, uninduced", xlab="Continuous Genome Position")
dev.off()
