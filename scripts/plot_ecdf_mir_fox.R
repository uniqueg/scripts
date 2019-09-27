## Alexander Kanitz
## 30-OCT-2014

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set working directory
setwd(args[1])

# Load fold change table (resulting from differential gene expression analysis)
de <- read.table("fold_changes.tab", stringsAsFactors=FALSE, skip=1, col.names=c("refseq_id", "logFC", "logCPM", "PValue"))

# Load Fox targets
fox <- read.table("fox_motifs_per_CLIPZ_transcript.tab", fill=TRUE, stringsAsFactors=FALSE, col.names=c("fox_motifs", "refseq_id"))

# Keep maximum number of motifs for duplicate entries
fox <- aggregate(fox_motifs~refseq_id, fox, max)

# Subset only expressed Fox targets
fox <- fox[fox$refseq_id %in% de$refseq_id, ]

# Add column indicating whether Fox target
fox <- cbind(fox, fox_target=fox$fox_motifs > 0)

# Merge fold changes and Fox targets
de_fox <- merge(de, fox)

# Load CLIPZ transcript ID <> Entrez gene ID conversion table
clipz_trx_to_gene <- read.table("clipz_trx_entrez_id_lookup_table", stringsAsFactors=FALSE, col.names=c("refseq_id", "entrez"))

# Add Entrez IDs by conversion of RefSeq IDs
last_col <- ncol(de_fox) + ncol(clipz_trx_to_gene) - 1
de_fox <- merge(de_fox, clipz_trx_to_gene)[ , c(last_col, 1:(last_col - 1))]

# Load miR targets file; skip unwanted header lines
mir_targets_hgnc_df <- read.csv("targets_of_fox_processed_miRs.csv", header=TRUE, sep=";", stringsAsFactors=FALSE, na.strings=c("NA", ""), skip=2)[-1,]

# Load Entrez gene ID <> HGNC symbol conversion table
entrez_2_hgnc <- read.delim("entrez_id_hgnc_symbol_lookup_table", stringsAsFactors=FALSE, col.names=c("entrez", "hgnc"))

# Convert miR target IDs to Entrez
mir_targets_entrez_ls <- sapply(mir_targets_hgnc_df, function(mir) sort(unique(subset(entrez_2_hgnc, subset=entrez_2_hgnc$hgnc %in% mir, select="entrez", drop=TRUE))))

# Indicate for each expressed transcript whether it is predicted to be targeted by each miR
de_fox_mir <- cbind(de_fox, lapply(mir_targets_entrez_ls, function(mir) de_fox$entrez %in% mir))

# Create pretty version for printing
de_fox_mir_print <- de_fox_mir
names(de_fox_mir_print) <- c("Entrez_gene_ID", "RefSeq_transcript_ID", "Log2_fold_change", "Log2_CPM", "P_value", "RBFOX2_motifs", "RBFOX2_target", gsub("\\.", "-", paste(names(de_fox_mir)[-1:-7], "target", sep="_")))

# Write out table
write.table(de_fox_mir_print, "entrez_refseq_fold_changes_p_values_fox_mir_targets", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


################### MIR ####################

# Group miR strands together
mir_pairs <- list(miR_20b=8:9, miR_21=10:11, miR_30=12:13, miR_32=14:15, miR_92=c(16,14), miR_106b=17, miR_181b=18:19, miR_107=20, miR_1_2=21)

# Create Fox vs strand target groups for each miR
groups <- lapply(mir_pairs, function(mir) {
        # Create empty list to hold results
        results <- list()
        # Check if there is more than one strand for current miRNA
        if (length(mir) == 1) {
                # From Fox target/non-target subsets, subset miR targets and non-targets
                results$mir_0 <- de_fox_mir[! de_fox_mir[, mir], 3]
                results$mir_1 <- de_fox_mir[de_fox_mir[, mir], 3]
        } else if (length(mir) == 2) {
                # From Fox target/non-target subsets, subset miR 5p, 3p and non-targets
                results$de_fox_mir_mir_5 <- de_fox_mir[de_fox_mir[, mir[1]], 3]
                results$de_fox_mir_mir_3 <- de_fox_mir[de_fox_mir[, mir[2]], 3]
                results$de_fox_mir_mir_0 <- de_fox_mir[! de_fox_mir$entrez %in% c(de_fox_mir[de_fox_mir[, mir[1]], 1], de_fox_mir[de_fox_mir[, mir[2]], 1]), 3]
                results$de_fox_mir_mir_5 <- de_fox_mir[de_fox_mir[, mir[1]], 3]
                results$de_fox_mir_mir_3 <- de_fox_mir[de_fox_mir[, mir[2]], 3]
                results$de_fox_mir_mir_0 <- de_fox_mir[! de_fox_mir$entrez %in% c(de_fox_mir[de_fox_mir[, mir[1]], 1], de_fox_mir[de_fox_mir[, mir[2]], 1]), 3]
        } else {
                # Debug
                print("THIS SHOULD NOT HAVE HAPPENED!")
        }
        # Return results
        return(results)
})
# Set list names
names(groups) <- names(mir_pairs)

# Get ECDF
ecdf <- lapply(groups, lapply, ecdf)

# Legends for plots
legend_6 <- c("5p targets w/o GCAUG", "3p targets w/o GCAUG", "Non-targets w/o GCAUG", "5p targets w/ GCAUG", "3p targets w/ GCAUG", "Non-targets w/ GCAUG")
legend_4 <- c("Non-targets w/o GCAUG", "Targets w/o GCAUG", "Non-targets w/ GCAUG", "Targets w/ GCAUG")




#################### FOX ###################

# Subset Fox targets by number of motifs
groups_fox <- list()
groups_fox$zero <- subset(de_fox_mir, subset=fox_motifs == 0, select=logFC, drop=TRUE)
groups_fox$one <- subset(de_fox_mir, subset=fox_motifs == 1, select=logFC, drop=TRUE)
groups_fox$two <- subset(de_fox_mir, subset=fox_motifs == 2, select=logFC, drop=TRUE)
groups_fox$three_to_four <- subset(de_fox_mir, subset=(fox_motifs >= 3 & fox_motifs <=4), select=logFC, drop=TRUE)
groups_fox$five_plus <- subset(de_fox_mir, subset=fox_motifs >= 5, select=logFC, drop=TRUE)
# Group names for output
pretty_group_names_fox <- c("0", "1", "2", "3 or 4", "5+")

# Get ECDF
ecdf_fox <- lapply(groups_fox, ecdf)

## Plot PDF
# Open graphics device
pdf("ecdf_by_number_of_rbfox2_motifs.pdf")
## Plot
plot(ecdf_fox$zero, ylim=c(0,1), xlim=c(-4,4), main="Fold change cumulative fractions by RBFOX2 motif count", xlab="Fold change siFox2 / siMock (log2)", ylab="Cumulative fraction", col=1, do.points=FALSE)
lines(ecdf_fox$one, col=2, do.points=FALSE)
lines(ecdf_fox$two, col=3, do.points=FALSE)
lines(ecdf_fox$three_to_four, col=4, do.points=FALSE)
lines(ecdf_fox$five_plus, col=5, do.points=FALSE)
legend(par("usr")[1], 1, legend=pretty_group_names_fox, col=1:5, lty=1, bty="n")
# Close graphics device
dump <- dev.off()

## Plot SVG
# Open graphics device
svg("ecdf_by_number_of_rbfox2_motifs.svg")
## Plot
plot(ecdf_fox$zero, ylim=c(0,1), xlim=c(-4,4), main="Fold change cumulative fractions by RBFOX2 motif count", xlab="Fold change siFox2 / siMock (log2)", ylab="Cumulative fraction", col=1, do.points=FALSE)
lines(ecdf_fox$one, col=2, do.points=FALSE)
lines(ecdf_fox$two, col=3, do.points=FALSE)
lines(ecdf_fox$three_to_four, col=4, do.points=FALSE)
lines(ecdf_fox$five_plus, col=5, do.points=FALSE)
legend(par("usr")[1], 1, legend=pretty_group_names_fox, col=1:5, lty=1, bty="n")
# Close graphics device
dump <- dev.off()


## Do two-sample Kolmogorov-Smirnov test
ks_mt <- suppressWarnings(sapply(groups_fox, function(group1) {
        ls <- list()
        for (group2_name in names(groups_fox)) {
                ls[[group2_name]] = ks.test(group1, groups_fox[[group2_name]])$p.value
        }
        return(ls)
}))
# Replace matrix dimnames with pretty names for output
dimnames(ks_mt) <- list(pretty_group_names_fox, pretty_group_names_fox)

# Write P value matrix to file
write.table(ks_mt, "stats.ks_test.fox_motifs.tab", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Get mean & median log2 fold changes
means <- sapply(groups_fox, mean)
medians <- sapply(groups_fox, median)

means_medians <- data.frame(motifs=pretty_group_names_fox, means=means, medians=medians)

write.table(means_medians, "stats.means_medians.fox_motifs.tab", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#################### MIR + FOX ####################

# Subset Fox targets and non-targets
fox_0 <- subset(de_fox_mir, subset=fox_target == FALSE)
fox_1 <- subset(de_fox_mir, subset=fox_target)

# Group miR strands together
mir_pairs <- list(miR_20b=8:9, miR_21=10:11, miR_30=12:13, miR_32=14:15, miR_92=c(16,14), miR_106b=17, miR_181b=18:19, miR_107=20, miR_1_2=21)

# Create Fox vs strand target groups for each miR
groups <- lapply(mir_pairs, function(mir) {
        # Create empty list to hold results
        results <- list()
        # Check if there is more than one strand for current miRNA
        if (length(mir) == 1) {
                # From Fox target/non-target subsets, subset miR targets and non-targets
                results$fox_0_mir_0 <- fox_0[! fox_0[, mir], 3]
                results$fox_0_mir_1 <- fox_0[fox_0[, mir], 3]
                results$fox_1_mir_0 <- fox_1[! fox_1[, mir], 3]
                results$fox_1_mir_1 <- fox_1[fox_1[, mir], 3]
        } else if (length(mir) == 2) {
                # From Fox target/non-target subsets, subset miR 5p, 3p and non-targets
                results$fox_0_mir_5 <- fox_0[fox_0[, mir[1]], 3]
                results$fox_0_mir_3 <- fox_0[fox_0[, mir[2]], 3]
                results$fox_0_mir_0 <- fox_0[! fox_0$entrez %in% c(fox_0[fox_0[, mir[1]], 1], fox_0[fox_0[, mir[2]], 1]), 3]
                results$fox_1_mir_5 <- fox_1[fox_1[, mir[1]], 3]
                results$fox_1_mir_3 <- fox_1[fox_1[, mir[2]], 3]
                results$fox_1_mir_0 <- fox_1[! fox_1$entrez %in% c(fox_1[fox_1[, mir[1]], 1], fox_1[fox_1[, mir[2]], 1]), 3]
        } else {
                # Debug
                print("THIS SHOULD NOT HAVE HAPPENED!")
        }
        # Return results
        return(results)
})
# Set list names
names(groups) <- names(mir_pairs)

# Get ECDF
ecdf <- lapply(groups, lapply, ecdf)

# Legends for plots
legend_6 <- c("5p targets w/o GCAUG", "3p targets w/o GCAUG", "Non-targets w/o GCAUG", "5p targets w/ GCAUG", "3p targets w/ GCAUG", "Non-targets w/ GCAUG")
legend_4 <- c("Non-targets w/o GCAUG", "Targets w/o GCAUG", "Non-targets w/ GCAUG", "Targets w/ GCAUG")

# SVG
invisible(lapply(names(ecdf), function(group_name) {

        # Re-format group, i.e. miRNA name
        mir_name <- gsub("\\_", "-", group_name)

        # Open plotting device
        svg(paste0("ecdf_", mir_name, ".svg"))

        # Plot eCDFs
        plot(ecdf[[group_name]][[1]], ylim=c(0,1), xlim=c(-4,4), main=mir_name, xlab="Fold change siFox2 / siMock (log2)", ylab="Cumulative fraction", do.points=FALSE, col=1)
        lines(ecdf[[group_name]][[2]], do.points=FALSE, col=2)
        lines(ecdf[[group_name]][[3]], do.points=FALSE, col=3)
        lines(ecdf[[group_name]][[4]], do.points=FALSE, col=4)

        # Check whether there are 6 or 4 curves; plot more lines and legend parameters accordingly
        if (length(ecdf[[group_name]]) == 6) {
                lines(ecdf[[group_name]][[5]], do.points=FALSE, col=5)
                lines(ecdf[[group_name]][[6]], do.points=FALSE, col=6)
                legend <- legend_6
                col <- 1:6
        } else if (length(ecdf[[group_name]]) == 4) {
                legend <- legend_4
                col <- 1:4
        } else {
                # Debug
                print("THIS SHOULD NOT HAVE HAPPENED!")
        }

        # Plot legend
        legend(par("usr")[1], 1, legend=legend, col=col, lty=1, bty="n")

        dev.off()

}))

# PDF
invisible(lapply(names(ecdf), function(group_name) {

        # Re-format group, i.e. miRNA name
        mir_name <- gsub("\\_", "-", group_name)

        # Open plotting device
        pdf(paste0("ecdf_", mir_name, ".pdf"))
 
        # Plot eCDFs
        plot(ecdf[[group_name]][[1]], ylim=c(0,1), xlim=c(-4,4), main=mir_name, xlab="Fold change siFox2 / siMock (log2)", ylab="Cumulative fraction", do.points=FALSE, col=1)
        lines(ecdf[[group_name]][[2]], do.points=FALSE, col=2)
        lines(ecdf[[group_name]][[3]], do.points=FALSE, col=3)
        lines(ecdf[[group_name]][[4]], do.points=FALSE, col=4)

        # Check whether there are 6 or 4 curves; plot more lines and legend parameters accordingly
        if (length(ecdf[[group_name]]) == 6) {
                lines(ecdf[[group_name]][[5]], do.points=FALSE, col=5)
                lines(ecdf[[group_name]][[6]], do.points=FALSE, col=6)
                legend <- legend_6
                col <- 1:6
        } else if (length(ecdf[[group_name]]) == 4) {
                legend <- legend_4
                col <- 1:4
        } else {
                # Debug
                print("THIS SHOULD NOT HAVE HAPPENED!")
        }

        # Plot legend
        legend(par("usr")[1], 1, legend=legend, col=col, lty=1, bty="n")

        dev.off()

}))

#### STATS

### MEANS
# Get mean log2 fold changes
means <- sapply(groups, lapply, mean)

## Convert mean and median "ragged" lists to filled matrices
nrow <- max(sapply(means, length))
means_mt <- sapply(means, function(mir) as.numeric(c(mir, rep(NA, nrow - length(mir)))))

# Add names for 6-member and 4-member groups
means_df <- data.frame(groups_6=legend_6, means_mt[,c(1,2,3,4,5,7,6,8,9)], groups_4=c(legend_4, rep(NA, 2)))

# Remove row names
rownames(means_df) <- NULL

# Re-format column names
colnames(means_df) <- c("2 miR strands", colnames(means_df)[c(-1,-length(colnames(means_df)))], "1 miR strand")

# Write out means and medians
write.table(means_df, "stats.means.tab", quote=FALSE, row.names=FALSE, col.names=TRUE, na="", sep="\t")

### MEDIANS
# Get median log2 fold changes
medians <- sapply(groups, lapply, median)

## Convert mean and median "ragged" lists to filled matrices
nrow <- max(sapply(medians, length))
medians_mt <- sapply(medians, function(mir) as.numeric(c(mir, rep(NA, nrow - length(mir)))))

# Add names for 6-member and 4-member groups
medians_df <- data.frame(groups_6=legend_6, medians_mt[,c(1,2,3,4,5,7,6,8,9)], groups_4=c(legend_4, rep(NA, 2)))

# Remove row names
rownames(medians_df) <- NULL

# Re-format column names
colnames(medians_df) <- c("2 miR strands", colnames(medians_df)[c(-1,-length(colnames(medians_df)))], "1 miR strand")

# Write out means and medians
write.table(medians_df, "stats.medians.tab", quote=FALSE, row.names=FALSE, col.names=TRUE, na="", sep="\t")

### KS-TEST
## Iterate over all groups/miRs
ks <- suppressWarnings(lapply(names(groups), function (group_name) {
        ## Do tests pairwise and return matrix of P values
        mt <- suppressWarnings(sapply(names(groups[[group_name]]), function(cond_1) {
                column <- suppressWarnings(sapply(names(groups[[group_name]]), function(cond_2) {
                        ks.test(groups[[group_name]][[cond_1]], groups[[group_name]][[cond_2]])$p.value
                }))
        }))
}))
# Set names for list elements/miRs
names(ks) <- names(groups)

## Pimp row and column names for printing
ks <- lapply(ks, function(mt) {
        if (nrow(mt) == 4) {
                row.names(mt) <- legend_4
                colnames(mt) <- legend_4
        } else if (nrow(mt) == 6) {
                row.names(mt) <- legend_6
                colnames(mt) <- legend_6
        }
        return(mt)
})

## Write out P value tables
for (group in names(ks)) {
        write.table(ks[[group]], paste("stats", "ks_test", group, "tab", sep="."), quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
}


