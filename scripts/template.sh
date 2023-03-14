#!/usr/bin/perl

#=============#
#  HEADER //  #
#=============#
## Created: ${date}
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: ${requirements}
#=============#
#  // HEADER  #
#=============#


#============#
#  USAGE //  #
#============#
usage()
### Returns u information for current script in a string
{
cat << USAGE
	Usage: $$0 ${short_usage_string}

	Description: ${script_description}

	Options:
	--${option}

	Comments:
	- ${usage_comment}

	Written by Alexander Kanitz, Biozentrum, Univerity of Basel on ${date}.
	Version 1.0 (${date})
USAGE
}
#-----------------------#
version ()
### Returns version information for current script in a string
{
cat <<VERSION;
$$0 version 1.0 (${date})
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on ${date}.
VERSION
}
#============#
#  // USAGE  #
#============#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES <---#
verbose=0
${other_options_variables}

#---> PARSE / ASSIGN OPTIONS <---#
while :
do
	case $$1 in
	--usage | --help)
		usage
		exit 0
		;;
	--version)
		version
		exit 0
		;;
	--verbose)
		verbose=1
		shift
		;;
	--) # End of all options
		shift
		break
		;;
	-*)
		echo "WARN: Unknown option (ignored): $$1" >&2
		shift
		;;
	*)  # no more options. Stop while loop
		break
		;;
	esac
done

#---> VERIFY OPTIONS <---#
# Print usage information and exit if required arguments are not specified
if [ ${required_argument} -eq '' ]; then 
	usage
	exit 0
fi

#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#
#  MAIN //  #
#===========#

#---> STATUS MESSAGE <---#
if [ $$verbose -eq 1 ]; then echo "Starting '$$0'...\n" >&2; fi

#---> MAIN VARIABLES <---#
${main_variables}

	#---> BODY <---#

	#---> ${main_code_block_comment} <---#
	${main_code_block}

#---> STATUS MESSAGE <---#
if [ $$verbose -eq 1 ]; then echo "Done.\n" >&2; fi

#---> PROGRAM EXIT <---#
exit 0;

#===========#
#  // MAIN  #
#===========#


#==================#
#  SUBROUTINES //  #
#==================#
sub ${sub_name} {
## Function: ${sub_function_description}
## Accepts: ${sub_argument_description}
## Returns: ${sub_return_value_description}
## Dependencies: ${sub_dependencies}
## Type: ${sub_type_generic_specific}

	#---> PASS ARGUMENTS ---#
	${passed_variables}

	#---> STATUS MESSAGE <---#
	if [ $$verbose -eq 1 ]; then echo "${sub_start_message}" >&2; fi

	#---> SUBROUTINE VARIABLES <---#
	${local_variables}

	#---> BODY <---#

		#---> ${sub_code_bloc} <---#
		${sub_code_block}

	#---> STATUS MESSAGE <---#
	if [ $$verbose -eq 1 ]; then echo "${sub_end_message}" >&2; fi

	#---> RETURN VALUE <---#
	return ${return_value};

}
#-----------------------#
${other_subs}
#==================#
#  // SUBROUTINES  #
#==================#