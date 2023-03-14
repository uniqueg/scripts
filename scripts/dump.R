##### CODE IN DEVELOPMENT

# FASTER METHOD FOR METADATA ADDITION
for (i in c(1:length(sc68_introns_list))) {
  print(i)
  if (length(sc68_introns_list[[i]]) > 0) {
	test <- sc68_introns_list[[i]]
    test[1]$source <- "TEST"#sc68_exons_list[[i]][1]$source
#    sc68_introns_list[[i]]$type <- "intron"
#    is.na(mcols(sc68_introns_list[[i]])) <- c("score","phase")
#    intron_no <- sc68_exons_list[[i]]$exon_number[1:(length(sc68_exons_list[[i]]$exon_number) - 1)]
#    sc68_introns_list[[i]]$exon_number <- paste0(intron_no, "/", intron_no + 1)
#    sc68_introns_list[[i]]$gene_biotype <- sc68_exons_list[[i]][1]$gene_biotype
#    sc68_introns_list[[i]]$gene_id <- sc68_exons_list[[i]][1]$gene_id
#    sc68_introns_list[[i]]$gene_name <- sc68_exons_list[[i]][1]$gene_name
#    is.na(mcols(sc68_introns_list[[i]])) <- "protein_id"
#    sc68_introns_list[[i]]$transcript_id <- sc68_exons_list[[i]][1]$transcript_id
#    sc68_introns_list[[i]]$transcript_name <- sc68_exons_list[[i]][1]$transcript_name    
#    print(paste0("INTRONS: ",length(sc68_introns_list[[i]])))
#    print(sc68_introns_list[[i]])

  }
}

# FASTER METHOD FOR METADATA ADDITION
# WORKS IN PRINCIPLE (<10 MIN), BUT: FINAL LIST CONTAINS ONLY LAST ELEMENT!!!
for (i in c(1:length(sc68_introns_list))) {
  print(i)
  if (length(sc68_introns_list[[i]]) > 0) {
    sc68_introns_list_test <- mapply(function(u,v) {	
      print(paste0("INTRON: ",length(v)))
      v$source <- u[1]$source
      v$type <- "intron"
      is.na(mcols(v)) <- c("score","phase")
      intron_no <- u$exon_number[1:(length(u$exon_number) - 1)]
      v$exon_number <- paste0(intron_no, "/", intron_no + 1)
      v$gene_biotype <- u[1]$gene_biotype
      v$gene_id <- u[1]$gene_id
      v$gene_name <- u[1]$gene_name
      is.na(mcols(v)) <- "protein_id"
      v$transcript_id <- u[1]$transcript_id
      v$transcript_name <- u[1]$transcript_name
      return(v)
    }, sc68_exons_list[i], sc68_introns_list[i])
  }
}

##### UNUSED CODE

# VALIDATION C1.
# Alternative method: use GenomicRanges::gaps method to infer
# 1. intergenic regions from "sc68_genes"
# 2. genic regions without overlaps from intergenic regions
sc68_intergenic_alt <- gaps(sc68_genes)
sc68_genic_alt <- gaps(sc68_intergenic_alt)
# Remove star (*) strand from "sc68_intergenic_alt"
sc68_intergenic_alt <- sc68_intergenic_alt[strand(sc68_intergenic_alt) != "*"]
# Differences in total interval lengths between "sc68_genic" and "sc68_genic_alt" should be 0
sum(width(sc68_genic)) - sum(width(sc68_genic_alt)) == 0
# TRUE!
# Differences in total interval lengths between "sc68_intergenic" and "sc68_intergenic_alt" should be 0
sum(width(sc68_intergenic)) - sum(width(sc68_intergenic_alt)) == 0
# TRUE!
# The sum of the total interval lengths of exonic and intronic should match the total interval length of genic regions
(sum(width(sc68_exonic)) + sum(width(sc68_intronic))) - sum(width(sc68_genic)) == 0
# TRUE!
# The sum of the total interval lengths of genic and intergenic regions should match the total length of the (double-stranded) genome
(sum(width(sc68_genic)) + sum(width(sc68_intergenic))) - (2 * sum(seqlengths(sc68_features))) == 0
# TRUE!
# Remove unused objects
rm(sc68_interexonic, sc68_genic_alt, sc68_intergenic_alt)

# VALIDATION C2.
# The sum of the total interval lengths of opposite exonic and opposite intronic should match the total interval length of opposite genic regions
(sum(width(sc68_opp_exonic)) + sum(width(sc68_opp_intronic))) - sum(width(sc68_opp_genic)) == 0
# TRUE!
# The sum of the total interval lengths of opposite exonic, opposite intronic and true intergenic should match the total interval length of intergenic regions
(sum(width(sc68_opp_exonic)) + sum(width(sc68_opp_intronic)) + sum(width(sc68_true_intergenic))) - sum(width(sc68_intergenic)) == 0
# TRUE!
# The sum of the total interval lengths of genic, opposite genic and true intergenic regions should match the total length of the (double-stranded) genome
(sum(width(sc68_genic)) + sum(width(sc68_opp_genic)) + sum(width(sc68_true_intergenic))) - (2 * sum(seqlengths(sc68_features))) == 0
# TRUE!

## GFF file export
# Remove names from "sc68_introns" to allow gff export
sc68_introns_red <- sc68_introns
names(sc68_introns_red) <- NULL
# NOT NECESSARY: convert intron numbers (1/2 to 1 etc)
mcols(sc68_introns_red)$exon_number <- as.integer(substr(mcols(sc68_introns_red)$exon_number,1,1))
# export files as gff
export(sc68_introns_v2, "sc68_introns_v2.gff", "gff")
export(sc68_introns, "sc68_introns.gff", "gff")
export(sc68_exons, "sc68_exons.gff", "gff")
export(sc68_genes, "sc68_genes.gff", "gff")
export(sc68_intergenic, "sc68_intergenic.gff", "gff")
# copy to local hdd (from local bash)
# scp akanitz@imlspenticton.uzh.ch:/home/akanitz/riboRNA/R_files/*.gff ~/Dropbox/UZH

# Return only those intervals that have an exon_number of > 1
sc68_multi_exons <- sc68_exons[sc68_exons$exon_number > 1, ]
# OPTIONAL: Check result 
sc68_multi_exons$exon_number
# Extract unique gene_id metadata entries as Character Vector
multi_exon_genes <- unique(mcols(sc68_multi_exons)$gene_id)
test[[2]][1] <- punion(test[[2]][2], test[[2]][1], fill.gap=TRUE)
index <- c(1:length(sc68_exons_list))
for (i in index)
{
	if (length(sc68_exons_list[[i]]) > 1)
		{	
		mini <- c(chrM=min(start(sc68_exons_list[[i]])))
		maxi <- c(chrM=max(end(sc68_exons_list[[i]])))
		intron <- gaps(sc68_exons_list[[i]], start=mini, end=maxi)		
		sc68_introns_list <- GRangesList(intron)		
		} 
}
intron <- gaps(sc68_exons_list[[2]], start=min(start(sc68_exons_list[[2]])), end=max(end(sc68_exons_list[[2]])))

sc68_introns_list_labels <- mendoapply(function(u,v,w) {
  print(length(w))
  if (length(w) > 0) {
    w$source <- u$source
    w$type <- "intron"
    is.na(mcols(w)) <- c("score","phase")
    intron_no <- v$exon_number[1:(length(v$exon_number) - 1)]
    w$exon_number <- paste0(intron_no, "/", intron_no + 1)
    w$gene_biotype <- u$gene_biotype
    w$gene_id <- u$gene_id
    w$gene_name <- u$gene_name
    is.na(mcols(w)) <- "protein_id"
    w$transcript_id <- u$transcript_id
    w$transcript_name <- u$transcript_name
  }
  return(w)
}, sc68_genes_list, sc68_exons_list, sc68_introns_list)

# Split sc68_exons to GenomicRangesList object > one list per gene (regardless of exon number)
sc68_exons_list <- split(sc68_exons, mcols(sc68_exons)$gene_id)
# !!!VERY SLOW!!!
# Define function that returns a GRanges object containing all genes as single entries (i.e. multiple exons are merged)
# 1. Generate a vector "index" with a numerical index for each "sc68_exons_list" entry
# 2. Loop through each "sc68_exons_list" entry via "index" (start for loop)
# 3. Assign first element of current "sc68_exons_list" entry to GRanges object "merged"
# 4. If multiple elements are present in current "sc68_exons_list" entry, compute the min/max start/end values of all elements and define as start/end of "merged"
# 5. If "index" is 1, assign "merged" to GRanges object "sc68_genes", else append it (end for loop)
# 6. Remove temporary variables
index <- c(1:length(sc68_exons_list))
for (i in index)
	{
	merged <- sc68_exons_list[[i]][1]
	if (length(sc68_exons_list[[i]]) > 1)
		{
		start(merged) <- min(start(sc68_exons_list[[i]]))
		end(merged) <- max(end(sc68_exons_list[[i]]))
		}
	if (i == 1)
		sc68_genes <- merged
	else	
		sc68_genes <- append(sc68_genes, merged)
	}
rm(i, index, merged)

## Infer introns
# Split sc68_exons to GenomicRangesList object > one list per gene (regardless of exon number)
sc68_exons_list <- split(sc68_exons, mcols(sc68_exons)$gene_id)
# !!!VERY SLOW!!!
# Define function that returns a GRanges object containing all introns
# 1. Generate a vector "index" with a numerical index for each "sc68_exons_list" entry
# 2. Loop through each "sc68_exons_list" entry via "index" (start outer for loop)
# 3. Check whether current list entry contains more than one element/interval (i.e. exon) (start outer if statement)
# 4. Generate a vector "index2" with a numerical index for each "sc68_exons_list[i]" GRanges object
# 5. Loop through each "sc68_exons_list[i]" entry via "index2" (start inner for loop)
# 6. Construct GRanges object "intron"; infer start/end values from flanking exons; derive all other data from first exon; set "score", "phase" and "protein_id" as NA
# 7. If "sc68_introns" is empty, assign "intron", else append it (start inner if statement, end inner/outer if statement, end inner/outer for loop)
# 8. Remove temporary variables
index <- c(1:length(sc68_exons_list))
for (i in index)
	{
	if (length(sc68_exons_list[[i]]) > 1)
		{
		index2 <- c(1:(length(sc68_exons_list[[i]])-1))
		for (j in index2)
			{
			intron <- GRanges(
			seqnames = seqnames(sc68_exons_list[[i]][1]),			
			ranges = IRanges(end(sc68_exons_list[[i]][j])+1,start(sc68_exons_list[[i]][j+1])-1),
			strand = strand(sc68_exons_list[[i]][1]),
			source = sc68_exons_list[[i]][j]$source,
			type = "intron",
			score = "",
			phase = "",
			exon_number = paste0(j,"/",j+1),
			gene_biotype = sc68_exons_list[[i]][1]$gene_biotype,
			gene_id = sc68_exons_list[[i]][1]$gene_id,
			gene_name = sc68_exons_list[[i]][1]$gene_name,
			protein_id = "",
			transcript_id = sc68_exons_list[[i]][1]$transcript_id,
			transcript_name = sc68_exons_list[[i]][1]$transcript_name
			)
			is.na(mcols(intron)) <- c("score","phase","protein_id")			
			if (i == 2 && j == 1)
				sc68_introns = as(intron, "GRanges")
			else 			
				sc68_introns = append(sc68_introns, intron) 			
			}		
		}
	}

#ug <- unlist(u)
#  
#setdiff(ug,vg)

rm(index, i, index2, j, intron)
# Set chromosome lengths (from http://www.yeastgenome.org/cache/genomeSnapshot.html # Genome Inventory table)
seqlengths(sc68_introns) <- c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779)
# Save
save(sc68_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "sc68_introns.rda"))
## VALIDATION / CROSS-COMPARISON WITH INTRONS
# Combine introns and exons
#sc68_exons_introns <- c(sc68_exons, sc68_introns)
# Remaining gaps should consitute intergenic regions
#sc68_intergenic_val <- gaps(sc68_exons_introns)
# PROBLEMS: 
# 3. Method is very slow!
# POSSIBLE SOLUTIONS:
# Overwrite/edit GenomicRanges::gaps function

# Use GenomicRanges::gaps function to fill in non-genic (i.e. intergenic) regions
sc68_intergenic <- gaps(sc68_genes)
# Add metadata
is.na(mcols(sc68_intergenic)) <- c("source","type","score","phase","exon_number","gene_biotype","gene_id","gene_name","protein_id","transcript_id","transcript_name")
sc68_intergenic$source <- "intergenic"
sc68_intergenic$type <- "intergenic"
sc68_intergenic$gene_biotype <- "intergenic"
# Save
save(sc68_intergenic, file=file.path("/home/akanitz/riboRNA/R_files/", "sc68_intergenic.rda"))
## VALIDATION / CROSS-COMPARISON WITH INTRONS
# Combine intergenic regions and exons
#sc68_exons_intergenic <- c(sc68_exons, sc68_intergenic)
# Remaining gaps should consitute introns
#sc68_introns_val <- gaps(sc68_exons_intergenic)
# PROBLEMS: 
# 1. Regions on the antisense strand of genes are regarded as intergenic (TRICKY)
# 2. Although not applicable/present in the yeast genome, the star strand (*) is considered (REMOVAL EASY)
# 3. Method is very slow!
# POSSIBLE SOLUTIONS:
# Overwrite/edit GenomicRanges::gaps function

## INFER INTERGENIC REGIONS FROM GENES AND EXONS
#
# Generate GRanges object with genes and their antisense strands
sc68_genes_ds <- c(sc68_genes, sc68_genes_as)
# Compute intergenic regions using the GRanges::gaps function
sc68_intergenic <- gaps(sc68_genes_ds)
# Remove star (*) strand from "sc68_intergenic" GRanges object
sc68_intergenic <- sc68_intergenic[strand(sc68_intergenic) != "*"]
# Add metadata
is.na(mcols(sc68_intergenic)) <- c("source","type","score","phase","exon_number","gene_biotype","gene_id","gene_name","protein_id","transcript_id","transcript_name")
sc68_intergenic$source <- "intergenic"
sc68_intergenic$type <- "intergenic"
sc68_intergenic$gene_biotype <- "intergenic"
# Save GRanges object "sc68_intergenic"
save(sc68_intergenic, file=file.path("/home/akanitz/riboRNA/R_files/", "sc68_intergenic.rda"))
##

# Obtain metadata from exons list >> Method 1: faster (<10 min), but less elegant
for (i in c(1:length(sc68_introns_list))) {
  print(i)
  if (length(sc68_introns_list[[i]]) > 0) {
    sc68_introns_list_test <- mapply(function(u,v) {	
      print(paste0("INTRON: ",length(v)))
      v$source <- u[1]$source
      v$type <- "intron"
      is.na(mcols(v)) <- c("score","phase")
      intron_no <- u$exon_number[1:(length(u$exon_number) - 1)]
      v$exon_number <- paste0(intron_no, "/", intron_no + 1)
      v$gene_biotype <- u[1]$gene_biotype
      v$gene_id <- u[1]$gene_id
      v$gene_name <- u[1]$gene_name
      is.na(mcols(v)) <- "protein_id"
      v$transcript_id <- u[1]$transcript_id
      v$transcript_name <- u[1]$transcript_name
      return(v)
    }, sc68_exons_list[i], sc68_introns_list[i])
  }
}

# ambiguous regions
sum(width(sc68_exons))
sc68_exons_excl <- setdiff(sc68_exons,sc68_introns)
sum(width(sc68_exons_excl))
sc68_exons_excl <- setdiff(sc68_exons_excl,sc68_exons_as)
sum(width(sc68_exons_excl))
sc68_exons_excl <- setdiff(sc68_exons_excl,sc68_introns_as)
sum(width(sc68_exons_excl))

sum(width(sc68_introns))
sc68_introns_excl <- setdiff(sc68_introns,sc68_exons)
sum(width(sc68_introns_excl))
sc68_introns_excl <- setdiff(sc68_introns_excl,sc68_exons_as)
sum(width(sc68_introns_excl))
sc68_introns_excl <- setdiff(sc68_introns_excl,sc68_introns_as)
sum(width(sc68_introns_excl))

sum(width(sc68_exons_as))
sc68_exons_as_excl <- setdiff(sc68_exons_as,sc68_exons)
sum(width(sc68_exons_as_excl))
sc68_exons_as_excl <- setdiff(sc68_exons_as_excl,sc68_introns)
sum(width(sc68_exons_as_excl))
sc68_exons_as_excl <- setdiff(sc68_exons_as_excl,sc68_introns_as)
sum(width(sc68_exons_as_excl))

sum(width(sc68_introns_as))
sc68_introns_as_excl <- setdiff(sc68_introns_as,sc68_exons)
sum(width(sc68_introns_as_excl))
sc68_introns_as_excl <- setdiff(sc68_introns_as_excl,sc68_introns)
sum(width(sc68_introns_as_excl))
sc68_introns_as_excl <- setdiff(sc68_introns_as_excl,sc68_exons_as)
sum(width(sc68_introns_as_excl))


sum(width(sc68_exons))
sum(width(sc68_introns))
sum(width(sc68_exons_as))
sum(width(sc68_introns_as))
print("===")
sum(width(sc68_exons_excl))
sum(width(sc68_introns_excl))
sum(width(sc68_exons_as_excl))
sum(width(sc68_introns_as_excl))

before <- c(sc68_exons, sc68_introns, sc68_exons_as, sc68_introns_as)
before <- gaps(before)
before <- gaps(before)
after <- c(sc68_exons_excl, sc68_introns_excl, sc68_exons_as_excl, sc68_introns_as_excl)
after <- gaps(after)
after <- gaps(after)

before <- c(sc68_exons, sc68_introns, sc68_exons_as, sc68_introns_as)	# all intervals
before <- gaps(before)							# gaps between intervals
before <- gaps(before)							# merged intervals


## C1. DEFINE GENIC (EXONIC + INTRONIC) AND INTERGENIC REGIONS
# Regions are exclusive, i.e. they do NOT contain overlaps; a priority scheme is applied to treat ambiguities
# Priority: Exonic > Intronic > Intergenic
# EXONIC
# To remove overlaps, use GenomicRanges::gaps method to infer
# 1. interexonic regions from "sc68_exons"
# 2. exonic regions without overlaps from interexonic regions
sc68_interexonic <- gaps(sc68_exons)
sc68_exonic <- gaps(sc68_interexonic)
save(sc68_exonic, file=file.path("~/riboRNA/R_files/", "sc68_exonic.rda"))
# EXCLUSIVE INTRONIC
# Use GenomicRanges::setdiff method to calculate the difference between "sc68_introns" and "sc68_exonic"
sc68_intronic <- setdiff(sc68_introns, sc68_exonic)
save(sc68_intronic, file=file.path("~/riboRNA/R_files/", "sc68_intronic.rda"))
# GENIC
# Append exonic and intronic regions
sc68_genic <- append(sc68_exonic, sc68_intronic)
save(sc68_genic, file=file.path("~/riboRNA/R_files/", "sc68_genic.rda"))
# INTERGENIC
# Use GenomicRanges::gaps method to infer intergenic from genic regions
sc68_intergenic <- gaps(sc68_genic)
# Remove star (*) strand from "sc68_intergenic"
sc68_intergenic <- sc68_intergenic[strand(sc68_intergenic) != "*"]
save(sc68_intergenic, file=file.path("~/riboRNA/R_files/", "sc68_intergenic.rda"))
##

## C2. ANNOTATE INTERGENIC REGIONS 
# Intergenic regions are split up into those regions that are or are not located opposite of exons, introns or genes (exons + introns) 
# Regions are exclusive, i.e. they do NOT contain overlaps; a priority scheme is applied to treat ambiguities
# Priority: opposite exonic > opposite intronic > "true intergenic"
# EXCLUSIVE OPPOSITE EXONIC
# Use GenomicRanges::setdiff to calculate the difference between "sc68_exonic_as" and "sc68_genic"
sc68_opp_exonic <- setdiff(sc68_exons_as, sc68_genic)
save(sc68_opp_exonic, file=file.path("~/riboRNA/R_files/", "sc68_opp_exonic.rda"))
# EXCLUSIVE OPPOSITE INTRONIC
# Use GenomicRanges::setdiff to calculate the difference between "sc68_exonic_as" and "sc68_genic"
sc68_opp_intronic <- setdiff(sc68_introns_as, append(sc68_genic, sc68_opp_exonic))
save(sc68_opp_intronic, file=file.path("~/riboRNA/R_files/", "sc68_opp_intronic.rda"))
# EXCLUSIVE OPPOSITE GENIC
# Append regions opposite exonic and intronic regions
sc68_opp_genic <- append(sc68_opp_exonic, sc68_opp_intronic)
save(sc68_opp_genic, file=file.path("~/riboRNA/R_files/", "sc68_opp_genic.rda"))
# "TRUE INTERGENIC"
# Use GenomicRanges::gaps method to infer true intergenic from genic and opposite genic regions
sc68_true_intergenic <- gaps(append(sc68_genic, sc68_opp_genic))
# Remove star (*) strand from "sc68_intergenic"
sc68_true_intergenic <- sc68_true_intergenic[strand(sc68_true_intergenic) != "*"]
save(sc68_true_intergenic, file=file.path("~/riboRNA/R_files/", "sc68_true_intergenic.rda"))
##


# GENIC
# Define exclusive genic (exonic, intronic and both) regions
sc68_genic <- c(sc68_excl_exonic, sc68_excl_intronic, sc68_exons_introns, ignore.mcols=TRUE)
save(sc68_genic, file=file.path("~/riboRNA/R_files/", "sc68_genic.rda"))
# INTERGENIC
# Use GenomicRanges::gaps method to infer intergenic from genic regions, remove star (*) strand from "sc68_intergenic" and save
sc68_intergenic <- gaps(sc68_genic)
sc68_intergenic <- sc68_intergenic[strand(sc68_intergenic) != "*"]
save(sc68_intergenic, file=file.path("~/riboRNA/R_files/", "sc68_intergenic.rda"))
##

# CDS
sc68_CDS <- sc68_features[values(sc68_features)[["type"]] == "CDS"]
save(sc68_CDS, file=file.path("~/riboRNA/R_files/", "sc68_CDS.rda"))
# START CODONS
sc68_start_codons <- sc68_features[values(sc68_features)[["type"]] == "start_codon"]
save(sc68_start_codons, file=file.path("~/riboRNA/R_files/", "sc68_start_codons.rda"))
# STOP CODONS
sc68_stop_codons <- sc68_features[values(sc68_features)[["type"]] == "stop_codon"]
save(sc68_stop_codons, file=file.path("~/riboRNA/R_files/", "sc68_stop_codons.rda"))

### E. CLEAN-UP
# Load all genome-related objects and save as image file
# Remove all objects
rm(list = ls())
# Create list of all saved genome-related objects
sc68_file_list = as.list(dir(pattern="^sc68.*rda$"))
# Load all elements in list
lapply(sc68_file_list, load, .GlobalEnv)
# Remove unused files
rm(sc68_file_list)
# Save image file of all genome-related objects
save.image(file = "sc68_objects_image.Rdata")
# Remove all objects
rm(list = ls())
###

# seqlengths(sc68_features) <- c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779)

# Remove all objects from environment
rm(list = ls())
# Load GenomicRanges package and genome-related objects
library(GenomicRanges)
load("sc68_objects_image.Rdata")
# Create list from all objects in environment
object_list <- as.list.environment(.GlobalEnv)
# Remove metadata columns
object_list <- lapply(object_list, function(u) {
  u <- u[, 0]
})
# Create GRangesList from genome-related objects
object_grl <- GRangesList(object_list)
# Add CUT/SUT intervals from Xu et al.
load("xu_objects_image.Rdata")
object_grl$xu_all_transcripts <- xu_all_transcripts
object_grl$xu_CUTs <- xu_CUTs
object_grl$xu_SUTs <- xu_SUTs
# Calculate total genome size (double-stranded)
genome_ds_size <- 2 * sum(seqlengths(sc68_features))
# Use sapply to acces object list
analysis_genome_features <- sapply(object_grl, function(u) {
  mt <- matrix()
  mt[1] <- sum(width(u))
  mt[2] <- round(mt[1, ] / genome_ds_size * 100,2)
  return(mt)
})
# Transpose results matrix
analysis_genome_features <- t(analysis_genome_features)
# Set matrix column names
colnames(analysis_genome_features) <- c("Total nt", "Total nt %")
# Sort by percentage of genome coverage (descending)
analysis_genome_features <- analysis_genome_features[order(analysis_genome_features[,"Total nt %"], decreasing = TRUE),]
# Save
save(analysis_genome_features, file=file.path("/home/akanitz/riboRNA/R_files/", "analysis_genome_features.Rdata"))
# Screen output
analysis_genome_features

file_names=as.list(dir(pattern="*.rda$"))
lapply(file_names,load, .GlobalEnv)
rm(file_names)
save.image(file = paste0(Sys.Date(),"_image.rda"))
library(GenomicRanges)

construct dataframe from this
sum(width(sc68_features))
sum(seqlengths(sc68_features))

> sum(seqlengths(sc68_features))
[1] 12157105
> sum(width(sc68_features))
[1] 18224768
> sum(width(sc68_CDS))
[1] 9030648
> sum(width(sc68_exons))
[1] 9153986
> sum(width(sc68_genes))
[1] 9270070
> sum(width(sc68_introns))
[1] 116084
> sum(width(sc68_intergenic))
[1] 27335892
> sum(width(sc68_ncRNA))
[1] 11438
> sum(width(sc68_protein_coding))
[1] 18121506
> sum(width(sc68_pseudo))
[1] 26367
> sum(width(sc68_rRNA))
[1] 27079
> sum(width(sc68_snoRNA))
[1] 13742
> sum(width(sc68_snRNA))
[1] 2408
> sum(width(sc68_start_codons))
[1] 20058
> sum(width(sc68_stop_codons))
[1] 20076
> sum(width(sc68_tRNA))
[1] 22228
> sum(width(xu_CUTs))
[1] 432590
> sum(width(xu_SUTs))
[1] 807567
> sum(width(xu_all_transcripts))
[1] 10917510
> sum(width(gerber_reads_A2))
[1] 1200262714
> sum(width(gerber_reads_A))
[1] 158658844
> sum(width(gerber_reads_B))
[1] 1136540267


## SAMPLE A1: Process bam
# convert bam to sam
samtools view -h S02_Ecogenics_120705_A-1-1.bam > riboA.sam;
# substitute chromosome entries in sam file
sed -i 's/S288c_chromosome/chr/' riboA.sam;
sed -i 's/S288c_mitochondrion/chrM/' riboA.sam;
# convert substituted sam to bam
samtools view -bS riboA.sam > riboA.bam;
# sort bam
samtools sort riboA.bam riboA_s;
# index bam (bai file)
samtools index riboA_s.bam;
# convert sorted, indexed bam to sam
samtools view -h riboA_s.bam > riboA_s.sam;

## SAMPLE A2: Process bam
# convert bam to sam
samtools view -h S02_Ecogenics_120705_A2-1-1.bam > riboA2.sam;
# substitute chromosome entries in sam file
sed -i 's/S288c_chromosome/chr/' riboA2.sam; 
sed -i 's/S288c_mitochondrion/chrM/' riboA2.sam;
# convert substituted sam to bam
samtools view -bS riboA2.sam > riboA2.bam;
# sort bam
samtools sort riboA2.bam riboA2_s;
# index bam (bai file)
samtools index riboA2_s.bam;
# convert sorted, indexed bam to sam
samtools view -h riboA2_s.bam > riboA2_s.sam;

## SAMPLE B: Process bam
# convert bam to sam
samtools view -h S02_Ecogenics_120705_B-1-1.bam > riboB.sam;
# substitute chromosome entries in sam file
sed -i 's/S288c_chromosome/chr/' riboB.sam; 
sed -i 's/S288c_mitochondrion/chrM/' riboB.sam;
# convert substituted sam to bam
samtools view -bS riboB.sam > riboB.bam;
# sort bam
samtools sort riboB.bam riboB_s;
# index bam (bai file)
samtools index riboB_s.bam;
# convert sorted, indexed bam to sam
samtools view -h riboB_s.bam > riboB_s.sam;

# GENIC
overlap_genic_reads <- summarizeOverlaps(sc68_genic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_genic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_genic_reads.rda"))
# INTERGENIC
overlap_intergenic_reads <- summarizeOverlaps(sc68_intergenic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_intergenic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_intergenic_reads.rda"))

######## NOT USED OR UNFINISHED


## CLEAN-UP
# Load all overlap objects and save as image file
# Remove all objects
rm(list = ls())
# Create list of all saved genome-related objects
overlap_file_list = as.list(dir(pattern="^overlap.*rda$"))
# Load all elements in list
overlap_object_list <- lapply(overlap_file_list, load, .GlobalEnv)
# Remove unused files
rm(overlap_file_list, overlap_object_list)
# Save image file of all genome-related objects
save.image(file = "overlap_objects_image.Rdata")
# Remove all objects
rm(list = ls())
###

### D. Clean-up
# > Load all genome-related objects and save as image files
# Remove all objects
rm(list = ls())
# Create list of all saved genome-related objects
sc68_file_list = as.list(dir(pattern="^sc68.*rda$"))
# Load all elements in list
sc68_object_list <- lapply(sc68_file_list, load, .GlobalEnv)
# Remove unused files
rm(sc68_file_list, sc68_object_list)
# Save image file of all genome-related objects
save.image(file = "sc68_objects_image.Rdata")
# Remove all objects
rm(list = ls())
###





sum(assays(exonic_reads_A_overlap)$counts)
sum(assays(intronic_reads_A_overlap)$counts)
sum(assays(opp_exonic_reads_A_overlap)$counts)
sum(assays(opp_intronic_reads_A_overlap)$counts)
sum(assays(true_intergenic_reads_A_overlap)$counts)

sum(assays(exonic_reads_A2_overlap)$counts)
sum(assays(intronic_reads_A2_overlap)$counts)
sum(assays(opp_exonic_reads_A2_overlap)$counts)
sum(assays(opp_intronic_reads_A2_overlap)$counts)
sum(assays(true_intergenic_reads_A2_overlap)$counts)

sum(assays(exonic_reads_B_overlap)$counts)
sum(assays(intronic_reads_B_overlap)$counts)
sum(assays(opp_exonic_reads_B_overlap)$counts)
sum(assays(opp_intronic_reads_B_overlap)$counts)
sum(assays(true_intergenic_reads_B_overlap)$counts)



## OVERLAP ANALYSIS
# GenomicRanges::summarizeOverlaps
# CDS
overlap_CDS_reads <- summarizeOverlaps(sc68_CDS, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_CDS_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_CDS_reads.rda"))
# EXONIC
overlap_exonic_reads <- summarizeOverlaps(sc68_exonic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_exonic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_exonic_reads.rda"))
# EXONS
overlap_exons_reads <- summarizeOverlaps(sc68_exons, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_exons_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_exons_reads.rda"))
# EXONS (ANTISENSE)
overlap_exons_as_reads <- summarizeOverlaps(sc68_exons_as, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_exons_as_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_exons_as_reads.rda"))
# FEATURES
overlap_features_reads <- summarizeOverlaps(sc68_features, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_features_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_features_reads.rda"))
# GENES
overlap_genes_reads <- summarizeOverlaps(sc68_genes, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_genes_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_genes_reads.rda"))
# GENES (ANTISENSE)
overlap_genes_as_reads <- summarizeOverlaps(sc68_genes_as, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_genes_as_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_genes_as_reads.rda"))
# GENIC
overlap_genic_reads <- summarizeOverlaps(sc68_genic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_genic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_genic_reads.rda"))
# INTERGENIC
overlap_intergenic_reads <- summarizeOverlaps(sc68_intergenic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_intergenic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_intergenic_reads.rda"))
# INTRONIC
overlap_intronic_reads <- summarizeOverlaps(sc68_intronic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_intronic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_intronic_reads.rda"))
# INTRONS
overlap_introns_reads <- summarizeOverlaps(sc68_introns, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_introns_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_introns_reads.rda"))
# INTRONS (ANTISENSE)
overlap_introns_as_reads <- summarizeOverlaps(sc68_introns_as, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_introns_as_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_introns_as_reads.rda"))
# NON-CODING RNA
overlap_ncRNA_reads <- summarizeOverlaps(sc68_ncRNA, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_ncRNA_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_ncRNA_reads.rda"))
# OPPOSITE OF EXONIC
overlap_opp_exonic_reads <- summarizeOverlaps(sc68_opp_exonic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_opp_exonic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_opp_exonic_reads.rda"))
# OPPOSITE OF GENIC
overlap_opp_genic_reads <- summarizeOverlaps(sc68_opp_genic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_opp_genic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_opp_genic_reads.rda"))
# OPPOSITE OF INTRONIC
overlap_opp_intronic_reads <- summarizeOverlaps(sc68_opp_intronic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_opp_intronic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_opp_intronic_reads.rda"))
# PROTEIN-CODING
overlap_protein_coding_reads <- summarizeOverlaps(sc68_protein_coding, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_protein_coding_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_protein_coding_reads.rda"))
# PSEUDOGENES
overlap_pseudo_reads <- summarizeOverlaps(sc68_pseudo, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_pseudo_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_pseudo_reads.rda"))
# RIBOSOMAL RNA
overlap_rRNA_reads <- summarizeOverlaps(sc68_rRNA, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_rRNA_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_rRNA_reads.rda"))
# SMALL NUCLEOLAR RNA
overlap_snoRNA_reads <- summarizeOverlaps(sc68_snoRNA, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_snoRNA_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_snoRNA_reads.rda"))
# SMALL NUCLEAR RNA
overlap_snRNA_reads <- summarizeOverlaps(sc68_snRNA, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_snRNA_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_snRNA_reads.rda"))
# START CODONS
overlap_start_codons_reads <- summarizeOverlaps(sc68_start_codons, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_start_codons_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_start_codons_reads.rda"))
# STOP CODONS
overlap_stop_codons_reads <- summarizeOverlaps(sc68_stop_codons, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_stop_codons_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_stop_codons_reads.rda"))
# TOTAL
overlap_total_reads <- summarizeOverlaps(sc68_total, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_total_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_total_reads.rda"))
# TRANSFER RNA
overlap_tRNA_reads <- summarizeOverlaps(sc68_tRNA, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_tRNA_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_tRNA_reads.rda"))
# TRUE INTERGENIC
overlap_true_intergenic_reads <- summarizeOverlaps(sc68_true_intergenic, list_bam_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_true_intergenic_reads, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_true_intergenic_reads.rda"))
##


## CLEAN-UP
# Load all overlap objects and save as image file
# Remove all objects
rm(list = ls())
# Create list of all saved genome-related objects
overlap_file_list = as.list(dir(pattern="^overlap.*rda$"))
# Load all elements in list
overlap_object_list <- lapply(overlap_file_list, load, .GlobalEnv)
# Remove unused files
rm(overlap_file_list, overlap_object_list)
# Save image file of all genome-related objects
save.image(file = "overlap_objects_image.Rdata")
# Remove all objects
rm(list = ls())
###

### D. Clean-up
# > Load all genome-related objects and save as image files
# Remove all objects
rm(list = ls())
# Create list of all saved genome-related objects
sc68_file_list = as.list(dir(pattern="^sc68.*rda$"))
# Load all elements in list
sc68_object_list <- lapply(sc68_file_list, load, .GlobalEnv)
# Remove unused files
rm(sc68_file_list, sc68_object_list)
# Save image file of all genome-related objects
save.image(file = "sc68_objects_image.Rdata")
# Remove all objects
rm(list = ls())
###

sum(assays(exonic_reads_A_overlap)$counts)
sum(assays(intronic_reads_A_overlap)$counts)
sum(assays(opp_exonic_reads_A_overlap)$counts)
sum(assays(opp_intronic_reads_A_overlap)$counts)
sum(assays(true_intergenic_reads_A_overlap)$counts)

sum(assays(exonic_reads_A2_overlap)$counts)
sum(assays(intronic_reads_A2_overlap)$counts)
sum(assays(opp_exonic_reads_A2_overlap)$counts)
sum(assays(opp_intronic_reads_A2_overlap)$counts)
sum(assays(true_intergenic_reads_A2_overlap)$counts)

sum(assays(exonic_reads_B_overlap)$counts)
sum(assays(intronic_reads_B_overlap)$counts)
sum(assays(opp_exonic_reads_B_overlap)$counts)
sum(assays(opp_intronic_reads_B_overlap)$counts)
sum(assays(true_intergenic_reads_B_overlap)$counts)

save(overlap_int_not_empty_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_all.rda"))
save(overlap_int_not_empty_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_CUTs.rda"))
save(overlap_int_not_empty_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_SUTs.rda"))
save(overlap_int_strict_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_SUTs.rda"))
save(overlap_int_strict_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_CUTs.rda"))
save(overlap_int_strict_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_all.rda"))
save(overlap_union_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_all.rda"))
save(overlap_union_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_SUTs.rda"))
save(overlap_union_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_CUTs.rda"))

load("overlap_int_not_empty_reads_xu_all.rda")
load("overlap_int_not_empty_reads_xu_CUTs.rda")
load("overlap_int_not_empty_reads_xu_SUTs.rda")
load("overlap_int_strict_reads_xu_SUTs.rda")
load("overlap_int_strict_reads_xu_CUTs.rda")
load("overlap_int_strict_reads_xu_all.rda")
load("overlap_union_reads_xu_all.rda")
load("overlap_union_reads_xu_SUTs.rda")
load("overlap_union_reads_xu_CUTs.rda")

## E1. LIST OF OBJECTS WITH OVERLAP TYPE "UNION"
# Remove all objects
rm(list = ls())
# Load elements
overlap_object_list <- lapply(list.files(pattern="^overlap_union.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_union_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_union_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_union_objects.rda"))
##

## E2. LIST OF OBJECTS WITH OVERLAP TYPE "INTERSECTION STRICT"
# Remove all objects
rm(list = ls())
# Load elements
overlap_object_list <- lapply(list.files(pattern="^overlap_int_strict.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_strict_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_strict_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_strict_objects.rda"))
##

## E3. LIST OF OBJECTS WITH OVERLAP TYPE "INTERSECTION NOT EMPTY"
# Remove all objects
rm(list = ls())
# Load elements
lapply(list.files(pattern="^overlap_int_not_empty.*rda$"), load, .GlobalEnv)
# Generate list from global environment
list_overlap_int_not_empty_objects <- as.list.environment(.GlobalEnv)
# Save list
save(list_overlap_int_not_empty_objects, file=file.path("~/riboRNA/R_files/", "list_overlap_int_not_empty_objects.rda"))
##


list <- as.list.environment(.GlobalEnv)

## C3. OTHER SC68 FEATURES
# Load objects
load("list_overlap_other_objects.rda")
# Use lapply to !!!!
other_reads <- lapply(seq_along(list_overlap_other_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_other_objects, n = gsub('overlap_int_not_empty_reads_', '', names(list_overlap_other_objects)))
 


## C3. OTHER SC68 FEATURES
# Load objects
load("list_overlap_other_objects.rda")
# Use lapply to !!!!
other_reads <- lapply(seq_along(list_overlap_other_objects), function(u, n, i) {
  sums <- matrix(ncol=3)
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
  colnames(sums) <- colnames(u[[i]])
  rownames(sums) <- n[[i]]
  return(sums)
}, u = list_overlap_other_objects, n = gsub('overlap_(union|int_strict|int_not_empty)_reads_', '', names(list_overlap_other_objects)))


union <- integer
int_strict <- integer
int_not_empty <- integer
for (i in seq_along(other_reads)) {
  if ('union_A2' %in% colnames(other_reads[[i]])[1]) union <- c(union, i)
  if ('int_strict_A2'  %in% colnames(other_reads[[i]])[1]) int_strict <- c(int_strict, i)
  if ('int_not_empty_A2'  %in% colnames(other_reads[[i]])[1]) int_not_empty <- c(int_not_empty, i)
}
other_reads_union <- list()
other_reads_union[[i]] <- other_reads[[i]]
other_reads_union <- other_reads[[union]]


do.call(cbind, lapply(other_reads, function (x) COMBINE(x, x, by = "row.names", all=TRUE))) 
   
     

lapply(other_reads, function(u) {
aggregate(u, by=colnames(u), FUN = "mean")
})


lapply(other_reads, function(u


do.call(rbind, test)

test <- lapply(other_reads, function(u) {
  x <- matrix()
  x <- COMBINE(x, u, by = "row.names", all=TRUE)
})
# Coerce to matrix
test2 <- as.matrix(test)





other_reads <- lapply(list_overlap_other_objects, function(u, v) {
  sums <- matrix(ncol=3)#(list(NULL, colnames(u)))
 # all <- matrix()#(list(NULL, colnames(u)))
  sums[,1] <- sum(assays(u[[i]][,1])$counts)
  sums[,2] <- sum(assays(u[[i]][,2])$counts)
  sums[,3] <- sum(assays(u[[i]][,3])$counts)
 print(colnames(u))
 print(v)
 #colnames(sums) <- v
 #rownames(sums) <- colnames(u)
 # all <- COMBINE(all, sums, by = "row.names", all=TRUE)
  return(sums)
}, v = gsub('overlap_(union|int_strict|int_not_empty)_reads_', '', names(list_overlap_other_objects)))#, simplify=FALSE)
# Coerce to matrix
other_reads <- do.call(cbind, other_reads)
other_reads

##
## C4. CUTS/SUTS
# Load objects
load("list_overlap_xu_objects.rda")
# Use lapply to obtain the total read counts for each genomic type in each sample
xu_reads <- lapply(list_overlap_xu_objects, function(u) assays(u)$counts)
# Coerce to matrix
xu_reads <- do.call(cbind, xu_reads)
##







## B3. COMBINE VALUES
genomic_analysis_mt <- cbind(genomic_type_annotation, genomic_type)
##
## B4. 


test <- sapply(analysis_genomic_type, function(u) {
  #print(u)
  cbind(u, u)
  #mt <- COMBINE(u, mt, by = "row.names", all = TRUE)
})




# Calculate total genome size (double-stranded)
genome_ds_size <- 2 * sum(seqlengths(list_sc68_genomic_type$exonic))
# Set matrix column names
colnames(analysis_genome_features) <- c("Total nt", "Total nt %")
# Sort by percentage of genome coverage (descending)
analysis_genome_features <- analysis_genome_features[order(analysis_genome_features[,"Total nt %"], decreasing = TRUE),]
# Save matrix
save(analysis_genome_features, file=file.path("/home/akanitz/riboRNA/R_files/", "analysis_genome_features.Rdata"))
# Screen output
analysis_genome_features
###



### C. GENE TYPE ANALYSIS








### C. FREQUENCIES IN GENOME
# Calculate the frequencies of different annotation types within the read files
## C1. OVERLAP MODE: "UNION"
# Load overlap objects
load("list_overlap_genomic_type_objects.rda")
# 
analysis_reads_union <- sapply(list_overlap_union_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[,1])$counts)
  mt[2] <- sum(assays(u[,2])$counts)
  mt[3] <- sum(assays(u[,3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_union <- t(analysis_reads_union)
# Set matrix column names
colnames(analysis_reads_union) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_union, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##
## C2. OVERLAP MODE: "INTERSECTION STRICT"
# Load overlap objects
load("list_overlap_int_strict_objects.rda")
# 
analysis_reads_int_strict <- sapply(list_overlap_int_strict_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_int_strict <- t(analysis_reads_int_strict)
# Set matrix column names
colnames(analysis_reads_int_strict) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_int_strict, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##
## C3. OVERLAP MODE: "INTERSECTION NOT EMPTY"
# Load overlap objects
load("list_overlap_int_not_empty_objects.rda")
# 
analysis_reads_int_not_empty <- sapply(list_overlap_int_not_empty_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_int_not_empty <- t(analysis_reads_int_not_empty)
# Set matrix column names
colnames(analysis_reads_int_not_empty) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_int_not_empty, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##






























# Remove all objects from environment
rm(list = ls())
# Load overlap objects
load("overlap_objects_image.Rdata")
# Create list from all objects in environment
overlap_object_list <- as.list.environment(.GlobalEnv)
# 
overlaps_genomic_analysis <- sapply(overlap_object_list, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
overlaps_genomic_analysis <- t(overlaps_genomic_analysis)
# Set matrix column names
colnames(overlaps_genomic_analysis) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save
save(overlaps_genomic_analysis, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))


# TO DO:
# Include relative numbers for read analysis
# Sort alphabetically BEFORE lapply (same order for both sc68 and overlap lists)
# Create character/factor vector for rownames
# Apply together with column names
# COMBINE matrices
# Sort whichever way seems best
# Save (genomic_analysis.Rdata)
# Plot (barplot)
# Statistics (chi square?)

list_features <- list_sc68_genomic_type
list_features$xu_all_transcripts <- xu_all_transcripts
list_features$xu_CUTs <- xu_CUTs
list_features$xu_SUTs <- xu_SUTs

load("overlap_int_not_empty_reads_exons.rda")
colnames(overlap_int_not_empty_reads_exons) <- paste("int_not_empty",basename(names(list_read_files)),sep="_")
colnames(overlap_int_not_empty_reads_exons) <- gsub('_s.bam','', colnames(overlap_int_not_empty_reads_exons))
colnames(overlap_int_not_empty_reads_exons) <- gsub('ribo','', colnames(overlap_int_not_empty_reads_exons))
save(overlap_int_not_empty_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_exons.rda"))


load("overlap_int_not_empty_reads_exons_as.rda")    
load("overlap_int_not_empty_reads_genes.rda")
load("overlap_int_not_empty_reads_genes_as.rda")   
load("overlap_int_not_empty_reads_gene_type.rda")
load("overlap_int_not_empty_reads_genomic_type.rda")
load("overlap_int_not_empty_reads_introns.rda")  
load("overlap_int_not_empty_reads_introns_as.rda")
load("overlap_int_not_empty_reads_xu_all.rda")
load("overlap_int_not_empty_reads_xu_CUTs.rda")
load("overlap_int_not_empty_reads_xu_SUTs.rda") 
load("overlap_int_strict_reads_exons.rda")
load("overlap_int_strict_reads_exons_as.rda")
load("overlap_int_strict_reads_genes.rda") 
load("overlap_int_strict_reads_genes_as.rda")
load("overlap_int_strict_reads_gene_type.rda")
load("overlap_int_strict_reads_genomic_type.rda")
load("overlap_int_strict_reads_introns.rda")
load("overlap_int_strict_reads_introns_as.rda")
load("overlap_int_strict_reads_xu_all.rda")
load("overlap_int_strict_reads_xu_CUTs.rda")
load("overlap_int_strict_reads_xu_SUTs.rda")

load("overlap_union_reads_exons.rda")
colnames(overlap_union_reads_exons) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_exons) <- gsub('_s.bam','', colnames(overlap_union_reads_exons))
colnames(overlap_union_reads_exons) <- gsub('ribo','', colnames(overlap_union_reads_exons))
save(overlap_union_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons.rda"))

load("overlap_union_reads_exons_as.rda")
colnames(overlap_union_reads_exons_as) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_exons_as) <- gsub('_s.bam','', colnames(overlap_union_reads_exons_as))
colnames(overlap_union_reads_exons_as) <- gsub('ribo','', colnames(overlap_union_reads_exons_as))
save(overlap_union_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons_as.rda"))

load("overlap_union_reads_genes.rda")
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))

load("overlap_union_reads_genomic_type.rda")
colnames(overlap_union_reads_genomic_type) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genomic_type) <- gsub('_s.bam','', colnames(overlap_union_reads_genomic_type))
colnames(overlap_union_reads_genomic_type) <- gsub('ribo','', colnames(overlap_union_reads_genomic_type))
save(overlap_union_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genomic_type.rda"))

load("overlap_union_reads_gene_type.rda")  
colnames(overlap_union_reads_gene_type) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_gene_type) <- gsub('_s.bam','', colnames(overlap_union_reads_gene_type))
colnames(overlap_union_reads_gene_type) <- gsub('ribo','', colnames(overlap_union_reads_gene_type))
save(overlap_union_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_gene_type.rda"))

load("overlap_union_reads_genes.rda")
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))

load("overlap_union_reads_introns.rda")
colnames(overlap_union_reads_introns) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_introns) <- gsub('_s.bam','', colnames(overlap_union_reads_introns))
colnames(overlap_union_reads_introns) <- gsub('ribo','', colnames(overlap_union_reads_introns))
save(overlap_union_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns.rda"))

load("overlap_union_reads_introns_as.rda")
colnames(overlap_union_reads_introns_as) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_introns_as) <- gsub('_s.bam','', colnames(overlap_union_reads_introns_as))
colnames(overlap_union_reads_introns_as) <- gsub('ribo','', colnames(overlap_union_reads_introns_as))
save(overlap_union_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns_as.rda"))

load("overlap_union_reads_xu_all.rda")
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))

load("overlap_union_reads_xu_CUTs.rda")
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))

load("overlap_union_reads_xu_SUTs.rda")
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))




overlap_union_reads_genes <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "Union", ignore.strand = FALSE)
colnames(overlap_union_reads_genes) <- paste("union",basename(names(list_read_files)),sep="_")
colnames(overlap_union_reads_genes) <- gsub('_s.bam','', colnames(overlap_union_reads_genes))
colnames(overlap_union_reads_genes) <- gsub('ribo','', colnames(overlap_union_reads_genes))



# GENE/EXON TYPES
overlap_union_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_gene_type.rda"))
# GENES
overlap_union_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_union_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_genes_as.rda"))
# EXONS
overlap_union_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_union_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_exons_as.rda"))
# INTRONS
overlap_union_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_union_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_introns_as.rda"))
##
## C2. MODE: INTERSECTION STRICT
# GENOMIC FEATURES
overlap_int_strict_reads_genomic_type <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genomic_type.rda"))
# GENE/EXON TYPES
overlap_int_strict_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_gene_type.rda"))
# GENES
overlap_int_strict_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_int_strict_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_genes_as.rda"))
# EXONS
overlap_int_strict_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_int_strict_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_exons_as.rda"))
# INTRONS
overlap_int_strict_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_int_strict_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_introns_as.rda"))
##
## C3. MODE: INTERSECTION NOT EMPTY
# GENOMIC FEATURES
overlap_int_not_empty_reads_genomic_type <- summarizeOverlaps(list_sc68_genomic_type, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_genomic_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genomic_type.rda"))
# GENE/EXON TYPES
overlap_int_not_empty_reads_gene_type <- summarizeOverlaps(list_sc68_gene_type, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_gene_type, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_gene_type.rda"))
# GENES
overlap_int_not_empty_reads_genes <- summarizeOverlaps(sc68_genes, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_genes, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genes.rda"))
# GENES (ANTISENSE)
overlap_int_not_empty_reads_genes_as <- summarizeOverlaps(sc68_genes_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_genes_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_genes_as.rda"))
# EXONS
overlap_int_not_empty_reads_exons <- summarizeOverlaps(sc68_exons, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_exons, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_exons.rda"))
# EXONS (ANTISENSE)
overlap_int_not_empty_reads_exons_as <- summarizeOverlaps(sc68_exons_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_exons_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_exons_as.rda"))
# INTRONS
overlap_int_not_empty_reads_introns <- summarizeOverlaps(sc68_introns, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_introns, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_introns.rda"))
# INTRONS (ANTISENSE)
overlap_int_not_empty_reads_introns_as <- summarizeOverlaps(sc68_introns_as, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_introns_as, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_introns_as.rda"))
###

### D. OVERLAP ANALYSIS CUTS/SUTS
# GenomicRanges::summarizeOverlaps
## D1. MODE: UNION
# CUTS
overlap_union_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_CUTs.rda"))
# SUTS
overlap_union_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_union_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "Union", ignore.strand = FALSE)
save(overlap_union_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_union_reads_xu_all.rda"))
##
## D2. MODE: INTERSECTION STRICT
# CUTS
overlap_int_strict_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_CUTs.rda"))
# SUTS
overlap_int_strict_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_int_strict_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "IntersectionStrict", ignore.strand = FALSE)
save(overlap_int_strict_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_strict_reads_xu_all.rda"))
##
## D3. MODE: INTERSECTION NOT EMPTY
# CUTS
overlap_int_not_empty_reads_xu_CUTs <- summarizeOverlaps(xu_CUTs, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_xu_CUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_CUTs.rda"))
# SUTS
overlap_int_not_empty_reads_xu_SUTs <- summarizeOverlaps(xu_SUTs, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_xu_SUTs, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_SUTs.rda"))
# ALL TRANSCRIPTS
overlap_int_not_empty_reads_xu_all <- summarizeOverlaps(xu_all_transcripts, list_read_files, mode = "IntersectionNotEmpty", ignore.strand = FALSE)
save(overlap_int_not_empty_reads_xu_all, file=file.path("/home/akanitz/riboRNA/R_files/", "overlap_int_not_empty_reads_xu_all.rda"))
###



require(grDevices) # for colours
tN <- table(Ni <- stats::rpois(100, lambda=5))
r <- barplot(tN, col=rainbow(20))
#- type = "h" plotting *is* 'bar'plot
lines(r, tN, type='h', col='red', lwd=2)

barplot(tN, space = 1.5, axisnames=FALSE,
        sub = "barplot(..., space= 1.5, axisnames = FALSE)")

barplot(VADeaths, plot = FALSE)
barplot(VADeaths, plot = FALSE, beside = TRUE)

mp <- barplot(VADeaths) # default
tot <- colMeans(VADeaths)
text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
barplot(VADeaths, beside = TRUE,
        col = c("lightblue", "mistyrose", "lightcyan",
                "lavender", "cornsilk"),
        legend = rownames(VADeaths), ylim = c(0, 100))
title(main = "Death Rates in Virginia", font.main = 4)

hh <- t(VADeaths)[, 5:1]
mybarcol <- "gray20"
mp <- barplot(hh, beside = TRUE,
        col = c("lightblue", "mistyrose",
                "lightcyan", "lavender"),
        legend = colnames(VADeaths), ylim= c(0,100),
        main = "Death Rates in Virginia", font.main = 4,
        sub = "Faked upper 2*sigma error bars", col.sub = mybarcol,
        cex.names = 1.5)
segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5)
stopifnot(dim(mp) == dim(hh))# corresponding matrices
mtext(side = 1, at = colMeans(mp), line = -2,
      text = paste("Mean", formatC(colMeans(hh))), col = "red")

# Bar shading example
barplot(VADeaths, angle = 15+10*1:5, density = 20, col = "black",
        legend = rownames(VADeaths))
title(main = list("Death Rates in Virginia", font = 4))

# border :
barplot(VADeaths, border = "dark blue") 


# log scales (not much sense here):
barplot(tN, col=heat.colors(12), log = "y")
barplot(tN, col=gray.colors(20), log = "xy")

# args.legend
barplot(height = cbind(x = c(465, 91) / 465 * 100,
                       y = c(840, 200) / 840 * 100,
                       z = c(37, 17) / 37 * 100),
        beside = FALSE,
        width = c(465, 840, 37),
        col = c(1, 2),
        legend.text = c("A", "B"),
        args.legend = list(x = "topleft"))
        
### B. LOAD FILES
# Genome files
load("list_sc68_genomic_type.rda")
# Remove metadata columns
list_sc68_genomic_type <- lapply(list_sc68_genomic_type, function(u) {
  u <- u[, 0]
})
# Gene files
load("list_sc68_gene_type.rda")
load("sc68_genes.rda")
load("sc68_genes_as.rda")
load("sc68_exons.rda")
load("sc68_exons_as.rda")
load("sc68_introns.rda")
load("sc68_introns_as.rda")
# Remove metadata columns
list_sc68_gene_type <- lapply(list_sc68_gene_type, function(u) {
  u <- u[, 0]
})
# CUT/SUT files; obtained from Xu et al. (see separate file "2_cut_sut_files.rda")
load("list_xu_objects.rda")
# Remove metadata columns
list_xu_objects <- lapply(list_xu_objects, function(u) {
  u <- u[, 0]
})
# Merge lists and add additional feature objects
list_features <- c(list_sc68_genomic_type, list_sc68_gene_type, list_xu_objects)
list_features$genes <- sc68_genes
list_features$genes_as <- sc68_genes_as
list_features$exons <- sc68_exons
list_features$exons_as <- sc68_exons_as
list_features$introns <- sc68_introns
list_features$introns_as <- sc68_introns_as
# Save
save(list_features, file=file.path("/home/akanitz/riboRNA/R_files/", "list_features.rda"))
###

### C. FREQUENCIES IN GENOME
# Calculate the frequencies of different annotation types in the genome
# Calculate total genome size (double-stranded)
load("sc68_features.rda")
genome_ds_size <- 2 * sum(seqlengths(sc68_features))
# Use sapply to access object list
analysis_genome_features <- sapply(list_features, function(u) {
  mt <- matrix()
  mt[1] <- sum(width(u))
  mt[2] <- round(mt[1, ] / genome_ds_size * 100,2)
  return(mt)
})
# Transpose results matrix
analysis_genome_features <- t(analysis_genome_features)
# Set matrix column names
colnames(analysis_genome_features) <- c("Total nt", "Total nt %")
# Sort by percentage of genome coverage (descending)
analysis_genome_features <- analysis_genome_features[order(analysis_genome_features[,"Total nt %"], decreasing = TRUE),]
# Save matrix
save(analysis_genome_features, file=file.path("/home/akanitz/riboRNA/R_files/", "analysis_genome_features.Rdata"))
# Screen output
analysis_genome_features
###

### C. FREQUENCIES IN GENOME
# Calculate the frequencies of different annotation types within the read files
## C1. OVERLAP MODE: "UNION"
# Load overlap objects
load("list_overlap_union_objects.rda")
# 
test <- names(list_overlap_union_objects)
analysis_reads_union <- sapply(list_overlap_union_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[,1])$counts)
  mt[2] <- sum(assays(u[,2])$counts)
  mt[3] <- sum(assays(u[,3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_union <- t(analysis_reads_union)
# Set matrix column names
colnames(analysis_reads_union) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_union, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##
## C2. OVERLAP MODE: "INTERSECTION STRICT"
# Load overlap objects
load("list_overlap_int_strict_objects.rda")
# 
analysis_reads_int_strict <- sapply(list_overlap_int_strict_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_int_strict <- t(analysis_reads_int_strict)
# Set matrix column names
colnames(analysis_reads_int_strict) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_int_strict, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##
## C3. OVERLAP MODE: "INTERSECTION NOT EMPTY"
# Load overlap objects
load("list_overlap_int_not_empty_objects.rda")
# 
analysis_reads_int_not_empty <- sapply(list_overlap_int_not_empty_objects, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
analysis_reads_int_not_empty <- t(analysis_reads_int_not_empty)
# Set matrix column names
colnames(analysis_reads_int_not_empty) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save matrix
save(analysis_reads_int_not_empty, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))
##

# Remove all objects from environment
rm(list = ls())
# Load overlap objects
load("overlap_objects_image.Rdata")
# Create list from all objects in environment
overlap_object_list <- as.list.environment(.GlobalEnv)
# 
overlaps_genomic_analysis <- sapply(overlap_object_list, function(u) {
  mt <- matrix()
  mt[1] <- sum(assays(u[, 1])$counts)
  mt[2] <- sum(assays(u[, 2])$counts)
  mt[3] <- sum(assays(u[, 3])$counts)
  return(mt)
})
# Transpose results matrix
overlaps_genomic_analysis <- t(overlaps_genomic_analysis)
# Set matrix column names
colnames(overlaps_genomic_analysis) <- c("A2 (transcriptome)", "A (translatome wt)", "B (translatome stress)")
# Save
save(overlaps_genomic_analysis, file=file.path("/home/akanitz/riboRNA/R_files/", "overlaps_genomic_analysis.rda"))

list_features <- list_sc68_genomic_type
list_features$xu_all_transcripts <- xu_all_transcripts
list_features$xu_CUTs <- xu_CUTs
list_features$xu_SUTs <- xu_SUTs
