# Set working directory
setwd("/import/bc2/home/zavolan/kanitz/PROJECTS/SpliceFactorsReprogramming")

# Global parameters
dataMatrix <- "analyzedData/kallistoQuant/kallistoQuant.tab"
groupAll <- c(1,rep(c(2,3,4),8))
outMDS500 <- "analyzedData/edgeR/MDS_500.pdf"
outMDS5000 <- "analyzedData/edgeR/MDS_5000.pdf"
sessionFile <- "analyzedData/edgeR/session.pdf"

# Load package
library(edgeR)

# Load data
df <- read.delim(dataMatrix, stringsAsFactors=FALSE)

# Functions
runEdgeR <- function(df, group, condName) {
    results <- list()
    dge_ls <- DGEList(df, group=group)
    dge_ls <- calcNormFactors(dge_ls)
    dge_ls <- estimateCommonDisp(dge_ls)
    dge_ls <- estimateTagwiseDisp(dge_ls)
    dge_ex <- exactTest(dge_ls)
    test_res <- decideTestsDGE(dge_ex)
    results$summ_tab <- summary(test_res)
    results$de <- sum(results$summ_tab[c(1,3)])
    return(results)
}
plot.MDS <- function(df, group, outFile, top=500) {
    dge_ls <- DGEList(df, group=group)
    pdf(file=outFile, width = 6, height = 6)
    plotMDS(dge_ls, top=top, labels=NULL, col=group)
    dev.off()
}

# Plot MDS
plot.MDS(df, groupAll, outMDS500, 500)
plot.MDS(df, groupAll, outMDS5000, 5000)

# Results container
results <- list()

#==========#

#---------------------------#
# ESRP1: Day 1-2 vs Day 7-8 #
#---------------------------#

# Local parameters
condName <- "ESRP_Days1-2_vs_Days_7-8"
colSelect <- c(3,6,21,24)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# Control: Day 1-2 vs Day 7-8 #
#-----------------------------#

# Local parameters
condName <- "CTRL_Days1-2_vs_Days_7-8"
colSelect <- c(2,5,20,23)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#---------------------------#
# DsRed: Day 1-2 vs Day 7-8 #
#---------------------------#

# Local parameters
condName <- "RED_Days1-2_vs_Days_7-8"
colSelect <- c(4,7,22,25)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#-----------------------------#
# ESRP1: Days 1-3 vs Days 6-8 #
#-----------------------------#

# Local parameters
condName <- "ESRP_Days1-3_vs_Days_6-8"
colSelect <- c(3,6,9,18,21,24)
group <- c(1,1,1,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-------------------------------#
# Control: Days 1-3 vs Days 6-8 #
#-------------------------------#

# Local parameters
condName <- "CTRL_Days1-3_vs_Days_6-8"
colSelect <- c(2,5,8,17,20,23)
group <- c(1,1,1,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# DsRed: Days 1-3 vs Days 6-8 #
#-----------------------------#

# Local parameters
condName <- "RED_Days1-3_vs_Days_6-8"
colSelect <- c(4,7,10,19,22,25)
group <- c(1,1,1,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#-----------------------------#
# ESRP1: Days 1-4 vs Days 5-8 #
#-----------------------------#

# Local parameters
condName <- "ESRP_Days1-4_vs_Days_5-8"
colSelect <- c(3,6,9,12,15,18,21,24)
group <- c(1,1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-------------------------------#
# Control: Days 1-4 vs Days 5-8 #
#-------------------------------#

# Local parameters
condName <- "CTRL_Days1-4_vs_Days_5-8"
colSelect <- c(2,5,8,11,14,17,20,23)
group <- c(1,1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# DsRed: Days 1-4 vs Days 5-8 #
#-----------------------------#

# Local parameters
condName <- "RED_Days1-4_vs_Days_5-8"
colSelect <- c(4,7,10,13,16,19,22,25)
group <- c(1,1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#-----------------------------#
# ESRP1: Days 0-2 vs Days 3-8 #
#-----------------------------#

# Local parameters
condName <- "ESRP_Days0-2_vs_Days_3-8"
colSelect <- c(1,3,6,9,12,15,18,21,24)
group <- c(1,1,1,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-------------------------------#
# Control: Days 0-2 vs Days 3-8 #
#-------------------------------#

# Local parameters
condName <- "CTRL_Days0-2_vs_Days_3-8"
colSelect <- c(1,2,5,8,11,14,17,20,23)
group <- c(1,1,1,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# DsRed: Days 0-2 vs Days 3-8 #
#-----------------------------#

# Local parameters
condName <- "RED_Days0-2_vs_Days_3-8"
colSelect <- c(1,4,7,10,13,16,19,22,25)
group <- c(1,1,1,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#-----------------------------#
# ESRP1: Days 0-2 vs Days 5-8 #
#-----------------------------#

# Local parameters
condName <- "ESRP_Days0-2_vs_Days_5-8"
colSelect <- c(1,3,6,15,18,21,24)
group <- c(1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-------------------------------#
# Control: Days 0-2 vs Days 5-8 #
#-------------------------------#

# Local parameters
condName <- "CTRL_Days0-2_vs_Days_5-8"
colSelect <- c(1,2,5,14,17,20,23)
group <- c(1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# DsRed: Days 0-2 vs Days 5-8 #
#-----------------------------#

# Local parameters
condName <- "RED_Days0-2_vs_Days_5-8"
colSelect <- c(1,4,7,16,19,22,25)
group <- c(1,1,1,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#-----------------------------#
# ESRP1: Days 0-1 vs Days 2-8 #
#-----------------------------#

# Local parameters
condName <- "ESRP_Days0-1_vs_Days_2-8"
colSelect <- c(1,3,6,9,12,15,18,21,24)
group <- c(1,1,2,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-------------------------------#
# Control: Days 0-1 vs Days 2-8 #
#-------------------------------#

# Local parameters
condName <- "CTRL_Days0-1_vs_Days_2-8"
colSelect <- c(1,2,5,8,11,14,17,20,23)
group <- c(1,1,2,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#-----------------------------#
# DsRed: Days 0-1 vs Days 2-8 #
#-----------------------------#

# Local parameters
condName <- "RED_Days0-1_vs_Days_2-8"
colSelect <- c(1,4,7,10,13,16,19,22,25)
group <- c(1,1,2,2,2,2,2,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#----------------------------#
# ESRP1 vs Control: Days 1-2 #
#----------------------------#

# Local parameters
condName <- "ESRP_vs_CTRL_Days1-2"
colSelect <- c(3,6,2,5)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#--------------------------#
# ESRP1 vs DsRed: Days 1-2 #
#--------------------------#

# Local parameters
condName <- "ESRP_vs_RED_Days1-2"
colSelect <- c(3,6,4,7)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#----------------------------#
# DsRed vs Control: Days 1-2 #
#----------------------------#

# Local parameters
condName <- "RED_vs_CTRL_Days1-2"
colSelect <- c(4,7,2,5)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#----------------------------#
# ESRP1 vs Control: Days 3-4 #
#----------------------------#

# Local parameters
condName <- "ESRP_vs_CTRL_Days3-4"
colSelect <- c(9,12,8,11)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#--------------------------#
# ESRP1 vs DsRed: Days 3-4 #
#--------------------------#

# Local parameters
condName <- "ESRP_vs_RED_Days3-4"
colSelect <- c(9,12,10,13)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#----------------------------#
# DsRed vs Control: Days 3-4 #
#----------------------------#

# Local parameters
condName <- "RED_vs_CTRL_Days3-4"
colSelect <- c(10,13,8,11)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#----------------------------#
# ESRP1 vs Control: Days 5-6 #
#----------------------------#

# Local parameters
condName <- "ESRP_vs_CTRL_Days5-6"
colSelect <- c(15,18,14,17)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#--------------------------#
# ESRP1 vs DsRed: Days 5-6 #
#--------------------------#

# Local parameters
condName <- "ESRP_vs_RED_Days5-6"
colSelect <- c(15,18,16,19)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#----------------------------#
# DsRed vs Control: Days 5-6 #
#----------------------------#

# Local parameters
condName <- "RED_vs_CTRL_Days5-6"
colSelect <- c(16,19,14,17)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#----------------------------#
# ESRP1 vs Control: Days 7-8 #
#----------------------------#

# Local parameters
condName <- "ESRP_vs_CTRL_Days7-8"
colSelect <- c(21,24,20,23)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#--------------------------#
# ESRP1 vs DsRed: Days 7-8 #
#--------------------------#

# Local parameters
condName <- "ESRP_vs_RED_Days7-8"
colSelect <- c(21,24,22,25)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#----------------------------#
# DsRed vs Control: Days 7-8 #
#----------------------------#

# Local parameters
condName <- "RED_vs_CTRL_Days7-8"
colSelect <- c(22,25,20,23)
group <- c(1,1,2,2)

# Subset data
df_select <- df[ , colSelect]

# Run edgeR
results[[condName]] <- runEdgeR(df_select, group, condName)

#==========#

#------------------------------------#
# Barplots: Number of DE transcripts #
#------------------------------------#
subset <- c("ESRP_Days1-2_vs_Days_7-8", "ESRP_Days1-3_vs_Days_6-8", "ESRP_Days1-4_vs_Days_5-8", "CTRL_Days1-2_vs_Days_7-8", "CTRL_Days1-3_vs_Days_6-8", "CTRL_Days1-4_vs_Days_5-8", "RED_Days1-2_vs_Days_7-8", "RED_Days1-3_vs_Days_6-8", "RED_Days1-4_vs_Days_5-8")
results_subset <- results[subset]
de_vec <- setNames(unlist(sapply(results_subset, "[", 2)), names(results_subset))
outFile <- "analyzedData/edgeR/start_vs_End.pdf"
pdf(file=outFile, width = 6, height = 6)
barplot(de_vec)
dev.off()
#----------------#
subset <- c("ESRP_Days0-1_vs_Days_2-8", "ESRP_Days0-2_vs_Days_3-8", "ESRP_Days0-2_vs_Days_5-8", "CTRL_Days0-1_vs_Days_2-8", "CTRL_Days0-2_vs_Days_3-8", "CTRL_Days0-2_vs_Days_5-8", "RED_Days0-1_vs_Days_2-8", "RED_Days0-2_vs_Days_3-8", "RED_Days0-2_vs_Days_5-8")
results_subset <- results[subset]
de_vec <- setNames(unlist(sapply(results_subset, "[", 2)), names(results_subset))
outFile <- "analyzedData/edgeR/start_vs_End_2.pdf"
pdf(file=outFile, width = 6, height = 6)
barplot(de_vec)
dev.off()
#----------------#
subset <- c("ESRP_vs_CTRL_Days1-2", "ESRP_vs_CTRL_Days3-4", "ESRP_vs_CTRL_Days5-6", "ESRP_vs_CTRL_Days7-8", "ESRP_vs_RED_Days1-2", "ESRP_vs_RED_Days3-4", "ESRP_vs_RED_Days5-6", "ESRP_vs_RED_Days7-8", "RED_vs_CTRL_Days1-2", "RED_vs_CTRL_Days3-4", "RED_vs_CTRL_Days5-6", "RED_vs_CTRL_Days7-8")
results_subset <- results[subset]
de_vec <- setNames(unlist(sapply(results_subset, "[", 2)), names(results_subset))
outFile <- "analyzedData/edgeR/conditions_vs_each_other.pdf"
pdf(file=outFile, width = 6, height = 6)
barplot(de_vec)
dev.off()

#--------------#
# Save session #
#--------------#
save.image(file=sessionFile)
