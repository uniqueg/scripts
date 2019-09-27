## Alexander Kanitz
## 18-DEC-2013

############ Pass #############
# 1. Coverage file            #
# 2. Prefix (e.g. library ID) #
###############################

# Get command-line arguments
args <- commandArgs(trailingOnly=TRUE)

## Define column names
col_names <- c("rseq", "source", "type", "start", "stop", "score", "str", "phase", "attr", "pos", "cov")

## Load data
df <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE, col.names=col_names)

## Extract names and IDs and parents
df$name <- sub("Name=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 3)))
df$id <- sub("ID=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 1)))
df$parent <- sub("Derives_from=", "", unlist(lapply(strsplit(df$attr, ";"), "[", 4)))

## Subset data by type
mat <- unique(df[df$type == "miRNA", c("start", "stop", "name", "id", "parent")])
pre <- unique(df[df$type == "miRNA_primary_transcript", c("rseq", "source", "type", "start", "stop", "score", "str", "phase", "attr", "name", "id")])
cov <- df[df$type == "miRNA_primary_transcript", c("name", "cov")]

## Convert precursor dataframe to list and add coverage
pre_ls <- lapply(split(pre, pre$name), as.list)
cov_ls <- split(cov, cov$name)
# Combine lists
for (name in names(cov_ls)) pre_ls[[name]]$cov <- cov_ls[[name]]$cov

## Add mature data to precursor list
pre_ls <- lapply(pre_ls, function(pre) {
	positions <- which(mat$parent == pre$id)
	pre$mature <- list()
	for (pos in positions) {
		pre$mature[[mat$name[pos]]] <- list(name=mat$name[pos], id=mat$id[pos], start=mat$start[pos], stop=mat$stop[pos])
	}
	return(pre)
})

## Plot
invisible(lapply(pre_ls, function(pre) { 
		
	# Get x coordinates
	x = seq(pre$start, pre$stop, 1)

	# Set labels
	main = paste0(args[2], ": ", pre$name, " / ", pre$id);
	xlab = pre$rseq
	ylab = "Coverage [nt]"

	# Open graphics device
	pdf(paste0(args[2], "_coverage_", pre$id, ".pdf"))
	
	# Plot coverage
	plot(x, pre$cov, type="l", main=main, xlab=xlab, ylab=ylab, xaxt="n")

	# Add custom x axis
	axis(1, at=c(pre$start, (pre$start + pre$stop) / 2, pre$stop))
	
	## Add mature miRNAs
	invisible(lapply(pre$mature, function(mat) {
		x_mat <- seq(mat$start, mat$stop, 1)
		y_mat <- rep(0, length(x_mat))
		lines(x_mat, y_mat, col="red", lwd=2, lty=2)
		text((mat$start + mat$stop) / 2, 0, labels=mat$name, pos=3, offset=0.5, col="red")
	}))
	
	# Close graphics device
	dev.off()

}))
