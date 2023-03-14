###############################################################################
#
# Author: herrmanch00
# Date:	14-MAY-2013
# Description: load, sort by FDR, and write sorted versions of edgeR output tables
#				"logFC_logCPM_P_FDR" from script edgeR_tables.R
###############################################################################

files <- dir(path="/home/herrmanch00/Documents/ALTSPLICE/checkGO/edgeR_output/tissues", pattern=glob2rx("inter*.tab"), full.names=TRUE, recursive=TRUE)

## Sort by FDR (ascending)
dump <- lapply(files, function(filename){
			t <- read.delim(filename)
			t <- t[order(t[,4]),]
			write.table(t,paste("sorted_by_FDR_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
			## FDR < 0.05
			s <- t[which(t[4] < 0.05),]
			write.table(s,paste("FDR_0.05_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
			## sort genes with FDR < 0.05 by logFC (descending)
			r <- s[order(-s[,1]),]
			write.table(r,paste("sorted_by_logFC_FDR_0.05_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
		})

## more stringent FDR (refseq_all_genes only)
files1 <- dir(path="/home/herrmanch00/Documents/ALTSPLICE/checkGO", pattern=glob2rx("inter*all*.tab"), full.names=TRUE, recursive=TRUE)

dump <- lapply(files1, function(filename){
			t <- read.delim(filename)
			t <- t[order(t[,4]),]
			## FDR < 1e-03
			s <- t[which(t[4] < 1e-03),]
			write.table(s,paste("FDR_1e-03_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
			## FDR < 1e-10
			s <- t[which(t[4] < 1e-10),]
			write.table(s,paste("FDR_1e-10_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
			## FDR top 1000
			s <- t[1:1000,]
			write.table(s,paste("FDR_top1000_",sub("inter_first_run_","",basename(filename)),sep=""), sep="\t")
			
		})