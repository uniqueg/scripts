# 30-JUL-2013
# Alexander Kanitz

usage="\nUsage: bash $0 [FASTA_FILE] [ID_FILE] [OUT_FILE] ([CHUNK_SIZE])\n\nDescription: Extracts entries (identifier and sequence line) with specified (partial) identifiers from a FASTA file\n\nRequired arguments:\n\n\t1. FASTA file; no line wrappings(!), i.e. one line per sequence, e.g.:\n\t\t>ENST00000012345\n\t\tACGTACGT\n\t\t>ENST00000011223\n\t\tGCGCATAT\n\t\t...\n\n\t2. TXT file of (partial) identifiers; one STRING per line, e.g.:\n\t\tENST00000012345\n\t\tENST00000011223\n\t\t...\n\n\t3. Output filename\n\nOptional arguments:\n\n\t1. Chunk size in lines [INT|DEFAULT: 1000]\n\nAttention: A temporary folder 'split' needs to be created in the working directory. The script will fail if there is already an existing folder of that name. The temporary files created in the folder as well as the folder itself will be deleted right before the script finishes, so make sure noone writes to it between creation and deletion!\n"

# Die if not all or too many arguments are provided
if [ "$#" -lt "3" -o "$#" -gt "4" ]
	then
		echo "[ERROR] Incorrect number of arguments."
		echo -e $usage
		exit
fi

# Die if FASTA file does not exist
if [ ! -f "$1" ]
	then
		echo "[ERROR] FASTA input file could not be opened."
		echo -e $usage
		exit
fi

# Die if ID file does not exist
if [ ! -f "$2" ]
        then
                echo "[ERROR] ID input file could not be opened."
                echo -e $usage
                exit
fi

# Set chunk size argument if provided, else use default
if [ "$4" ]
	then 
		# Die if specified chunk size is not an integer
		if ! echo $4 | egrep -q '^[0-9]+$'; then
			echo "[ERROR] Illegal chunk size. Use integer."
                	echo -e $usage
                	exit
		fi
		chunk=$4
	else
		chunk=1000
fi

# Die if folder './split' exists in working directory, else create...
if [ -d ./split ]
	then
		echo "[ERROR] Temporary directory could not be created." 
		echo -e $usage
		exit
	else
		mkdir ./split
fi

# Write status message
echo -e "Splitting identifier file '$2'..."
# Split ID file into chunks of 5000 lines
split -l $chunk -a 4 $2 ./split/ids_

# Write status message
echo "Extracting entries..."
# For each file chunk...
for file in ./split/ids_*
	do
		# Grep lines containing current chunk of identifiers as well as the following sequence lines and write them to a temporary output file
		grep -A1 --no-group-separator -f $file $1 > ${file}_seq
		# Write status message
		echo -e "\tWritten temporary FASTA file '${file}_seq'."
	done

# Write status message
echo "Merging temporary FASTA files to output FASTA file '$3'..."
# Merge temporary output files to output file
cat ./split/ids_*_seq > $3

# Write status message
echo "Cleaning up..."
# Remove temporary files in folder './split', then attempt to delete './split' itself
rm -f ./split/ids_*
rmdir ./split

# Write status message
echo "Done."
