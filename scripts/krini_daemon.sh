#!/import/bc2/home/zavolan/krini/bin/bash
# Written by Alexander Kanitz
# 17-SEP-2013
# v0.9

# DESCRIPTION: The daemon watches a specified folder for incoming Anduril run files, moves them to dedicated folders, executes them in the background and updates the specified database. Afterwards it watches the corresponding processes until they finish and updates the status field in the database accordingly.

# USAGE: bash ./krini_daemon.sh

# REQUIREMENTS: Bash v4.0+ and sqlite v2.1+ (but smaller than v3!) must be in $PATH

# NOTES:Verify/set parameters in the declaration section in PART A below


#---> PART A: PRE-REQUISITES <---#
	
	## Declarations
	startup="$HOME/.bashrc"				# Config file to execute at start
	watch_folder="$HOME/.WATCH"			# Folder to watch
	db_file="$HOME/KriniDatabase/KriniDB.db"	# SQLite database file
	interval=3					# Time (in sec) b/w checking processes
	status_started="PROCESSING"	# Status code for started Anduril run
	status_error="ERROR"		# Status code for Anduril run that finished w/ error(s)
	status_completed="DONE"		# Status code for Anduril run that finished w/o errors
	
	# Source startup file
	source $startup
	
	# Declare associative array
	declare -A started_processes
	
	## Set input field separator to line breaks and carriage returns only
	SAVEIFS=$IFS
	IFS=$(echo -en "\n\b")
	
#---> END PART A <---#



#---> PART B: START ANDURIL PROCESSES <---#

## Iterate over every file in folder $HOME/.WATCH..
	for file in $watch_folder/*
	do
		# Extract basename of current file
		base=`basename $file`
	
		# Skip if file is info file
		if [[ "$base" == "DO_NOT_PLACE_FILES_HERE" ]]
		then
			continue
	
		## If file looks like proper Anduril run command file (filename starts with "execute_", file contains a line starting "run_dir='" and a line starting "anduril run")...
		elif [[ "$base" == execute_* ]] #&& grep -q "^anduril_run" "$file" && grep -q "^run_dir='" "$file"
		then
			# Extract run directory (format: run_dir='DIR/run_serial_number') and run serial number
			run_dir=`cat "$file" | awk -F'=' '{if ($1 == "run_dir") {print $2}}' | awk -F\' '{print $2}'`
			run_sn=`basename $run_dir`

			# Move file to run directory and update $file variable with new path
			mv "$file" "$run_dir/LOG"
			file="${run_dir}/LOG/${base}"

			# Execute file in background
			bash "$file" &> /dev/null &

			# Load "process ID" -> "run serial number" pair into associative array
			started_processes["$!"]="$run_sn"

			## Update status for current run in SQLite database file
			update_db=`echo "UPDATE experiment SET status = '$status_started' WHERE ex_id = $run_sn"`
			sqlite $db_file "$update_db"
	
		## Otherwise remove file / directory
		else
			rm -r -f "$file"
		fi
	done
	
#---> END PART B <---#



#---> PART C: WATCH STARTED ANDURIL PROCESSES <---#
	
	## As long as one or more started processes are still running...
	while [ "${#started_processes[@]}" -ne 0 ]
	do
		## Sleep for specified interval
		sleep $interval
	
		## Check the status of each process ID...
		for pid in "${!started_processes[@]}"
		do
			## If process finished...
			if [[ `ps -p $pid -o pid=` -ne $pid ]]
			then
				## Detect process change and get exit status
				wait $pid > /dev/null
				exit_code=$?

				## Update status for current run serial number in SQLite database file
				if [[ $exit_code -eq 0 ]] 
				then
	                                update_db=`echo "UPDATE experiment SET status = '$status_completed' WHERE ex_id = ${started_processes["$pid"]};"`
	                                sqlite $db_file "$update_db"
				else
	                                error="$status_error: $exit_code"
	                                update_db=`echo "UPDATE experiment SET status = '$error' WHERE ex_id = ${started_processes["$pid"]};"`
	                                sqlite $db_file "$update_db"
				fi

				# Remove "process ID" -> "run serial number" pair
				unset started_processes["$pid"]
			fi
		done
	done
	
#---> END PART C <---#



#---> PART D: CLEANUP <---#
	
	# Reset input field separator to original value
	IFS=$SAVEIFS
	
	# Return with zero exit status
	exit 0
	
#---> END PART D  <---#
