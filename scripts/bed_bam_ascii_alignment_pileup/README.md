# ASCII-style alignment pileup

## Description

Generates an ASCII-style pileup of read alignments in one or more BAM files against one or 
more regions specified in a BED file.

## Help screen

```console
Usage: ./bed_bam_ascii_alignment_pileup.R [OPTIONS] --bed <PATH> --bam <PATH>

Generates an ASCII-style pileup alignment of read alignments in one or more BAM
files against one or more regions specified in a BED file.

Author: Alexander Kanitz
Affiliation: Biozentrum, University of Basel
Email: alexander.kanitz@alumni.ethz.ch
Version: 1.0.0 (28-SEP-2019)
Requires: optparse, rtracklayer, GenomicFeatures, GenomicAlignments

Options:
	--bed=FILE
		Required: Path to BED file of query regions. Consider the
            `--maximum-region-width` parameter.

	--bam=FILE
		Required: Path to BAM file containing the alignments to be visualized against
            the query region(s). It is possible to use glob wildcards to select multiple BAM files
            at once. However, the feature is experimental and may lead to inconsistencies if a
            read contains an indel.

	--reference=FILE
		Reference genome sequence in FASTA format. The file *MUST* be compressed with
            BGZIP. If supplied, the reference sequence for the query region(s) will be added to
            the output. Note that on the first run with a specific reference genome file, an FAI
            index is generated which will take some time.

	--annotations=FILE
		Annotation file in GFF/GTF format used to annotate sequences. If supplied,
            features overlapping the query region(s) will be visualized in the output. Ensure that
            the argument to option `annotation-name-field` corresponds to a field in the
            annotations, otherwise the script will fail.

	--output-directory=DIR
		Output directory. One output file will be created for each region in `--bed` and
            the filenames will be generated from the basenames of the supplied BAM file(s) and the
            name field (4th column) of the BED file. [default "."]

	--maximum-region-width=INT
		Maximum input region width. Use with care as wide regions will use excessive
            resources. [default 200]

	--do-not-collapse-alignments
		Show alignments of reads with identical sequences individually.

	--minimum-count=INT
		Alignments of reads with less copies than the specified number will not be
            printed. Option is not considered if `do-not-collapse-alignments` is set.
            [default 1]

	--annotation-name-field=STR
		Annotation field used to populate the `name` column in the output.
            [default "Name"]

	--padding-character=CHAR
		Character used for padding alignments. [default .]

	--indel-character=CHAR
		Character to denote insertions and deletions in alignments. [default -]

	-h, --help
		Show this information and die.

	-v, --verbose
		Print log messages to STDOUT.

```

## Example output

```console
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	hsa-let-7a-1
....>>>>>>>>>>>>>>>>>>>>>>.....................................................	hsa-let-7a-5p
.......................................................>>>>>>>>>>>>>>>>>>>>>...	hsa-let-7a-3p
GGGATGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTCCTA	9:94175958-94176036:+
..CATGAGGTAGTAGGTTGTATAGT......................................................	1
...A-GAGGTAGTAGGTTGTATAG.......................................................	410
...A-GAGGTAGTAGGTTGTATAGTT.....................................................	129
...A-GAGGTAGTAGGTTGTATAGT......................................................	106
...CTGAGGTAGTAGGTTGTATAG.......................................................	18
...CTGAGGTAGTAGGTTGTATAGT......................................................	3
...A-GAGGTAGTAGGTTGTA..........................................................	2
...A-GAGGTAGTAGGTTGTAT.........................................................	1
...AT-AGGTAGTAGGTTGTATAG.......................................................	1
...CTGAGGTAGTAGGTTGTATAGTT.....................................................	1
....TGAGGTAGTAGGTTGTATAG.......................................................	71019
....TGAGGTAGTAGGTTGTATAGTT.....................................................	22757
....TGAGGTAGTAGGTTGTATAGT......................................................	18148
....TGAGGTAGTAGGTTG............................................................	3367
....TGAGGTAGTAGGTTGTATAGTTT....................................................	700
....TGAGGTAGTAGGTTGTATA........................................................	380
(truncated)
```
