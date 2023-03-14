## PLOT YH2AX CLUSTERS AND ASISI SITES 

# Import .bed file as GRanges object
library(rtracklayer)
yH2Ax <- import("clusters_100000_z_cutoff_induced.bed", genome="hg19", asRangedData=FALSE)
AsiSI <- import("AsiSI_orig.bed", genome="hg19", asRangedData=FALSE) 

# Extract genome and chromosome information
gen <- as.character(unique(genome(c(yH2Ax,AsiSI))))
chr <- names(genome(yH2Ax))
chr <- chr[!grepl("_", chr)]		# Remove incomplete/unknown chromosomes
chr <- chr[!grepl("chrM", chr)]		# Remove mitochondrial chromosome

# Gviz
library(Gviz)
gTrack <- GenomeAxisTrack()
iTrack <- IdeogramTrack(genome=gen, chromosome=chr, showId=FALSE, fill=FALSE, col=FALSE)
yH2AxTrack <- AnnotationTrack(yH2Ax, name="yH2Ax clusters", stacking="dense")
AsiSITrack <- AnnotationTrack(AsiSI, name="AsiSI recognition sites", stacking="dense", fill="orange")

lapply(chr, function(x) {
	#pdf(paste(x,".pdf", sep=""), width=12, height=3)
	plotTracks(list(iTrack, gTrack, yH2AxTrack, AsiSITrack), from=NULL, to=NULL, chromosome=x)
	#dev.off()
})


## PLOT YH2AX CLUSTERS AND ASISI SITES 

# Import .bed file as GRanges object
library(rtracklayer)
yH2Ax_i <- import("clusters_100000_z_cutoff_induced.bed", genome="hg19", asRangedData=FALSE)
yH2Ax_u <- import("clusters_100000_z_cutoff_uninduced.bed", genome="hg19", asRangedData=FALSE)
AsiSI <- import("AsiSI_orig.bed", genome="hg19", asRangedData=FALSE)
yH2Ax_AsiSI_i <- import("overlap_clusters_100000_induced_AsiSI.bed", genome="hg19", asRangedData=FALSE)
yH2Ax_AsiSI_u <- import("overlap_clusters_100000_uninduced_AsiSI.bed", genome="hg19", asRangedData=FALSE)
Overlap_iu <- import("overlap_clusters_100000_induced_vs_uninduced.bed", genome="hg19", asRangedData=FALSE)
Iacovoni_i <- import("yH2Ax_domains_Iacovoni_EMBO_hg19_induced.bed", genome="hg19", asRangedData=FALSE)
Iacovoni_u <- import("yH2Ax_domains_Iacovoni_EMBO_hg19_uninduced.bed", genome="hg19", asRangedData=FALSE)
Overlap_Iacovoni_i <- import("overlap_ChIP_with_Iacovoni_EMBO_induced.bed", genome="hg19", asRangedData=FALSE)
Overlap_Iacovoni_u <- import("overlap_ChIP_with_Iacovoni_EMBO_uninduced.bed", genome="hg19", asRangedData=FALSE)

# Extract genome and chromosome information
gen <- as.character(unique(genome(AsiSI)))
chr <- names(genome(AsiSI))
chr <- chr[!grepl("_", chr)]		# Remove incomplete/unknown chromosomes
chr <- chr[!grepl("chrM", chr)]		# Remove mitochondrial chromosome

# Gviz
library(Gviz)
gTrack <- GenomeAxisTrack()
iTrack <- IdeogramTrack(genome=gen, chromosome=chr, showId=FALSE, fill=FALSE, col=FALSE)
yH2Ax_i_Track <- AnnotationTrack(yH2Ax_i, name="yH2Ax clusters, induced", stacking="dense")
yH2Ax_u_Track <- AnnotationTrack(yH2Ax_u, name="yH2Ax clusters, uninduced", stacking="dense", fill="orange")
AsiSI_Track <- AnnotationTrack(AsiSI, name="AsiSI recognition sites", stacking="dense", fill="green")
yH2Ax_AsiSI_i_Track <- AnnotationTrack(yH2Ax_AsiSI_i, name="yH2Ax clusters with AsiSI sites, induced", stacking="dense")
yH2Ax_AsiSI_u_Track <- AnnotationTrack(yH2Ax_AsiSI_u, name="yH2Ax clusters with AsiSI sites, uninduced", stacking="dense", fill="orange")
Overlap_iu_Track <- AnnotationTrack(Overlap_iu, name="Overlap induced vs. uninduced", stacking="dense", fill="red")
Iacovoni_i_Track <- AnnotationTrack(Iacovoni_i, name="Iacovoni et al. clusters, induced", stacking="dense")
Iacovoni_u_Track <- AnnotationTrack(Iacovoni_u, name="Iacovoni et al. clusters, uninduced", stacking="dense", fill="orange")
Overlap_Iacovoni_i_Track <- AnnotationTrack(Overlap_Iacovoni_i, name="Overlap with Iacovoni et al., induced", stacking="dense", fill="red")
Overlap_Iacovoni_u_Track <- AnnotationTrack(Overlap_Iacovoni_u, name="Overlap with Iacovoni et al., uninduced", stacking="dense", fill="red")

lapply(chr, function(x) {
	pdf(paste(x,".pdf", sep=""), width=12, height=6)
	plotTracks(list(iTrack, gTrack, yH2Ax_i_Track, yH2Ax_u_Track, AsiSI_Track, yH2Ax_AsiSI_i_Track, yH2Ax_AsiSI_u_Track, Overlap_iu_Track, Iacovoni_i_Track, Iacovoni_u_Track, Overlap_Iacovoni_i_Track, Overlap_Iacovoni_u_Track), from=NULL, to=NULL, chromosome=x)
	dev.off()
})