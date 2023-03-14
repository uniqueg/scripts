#!/usr/bin/perl

#=============#
#  HEADER //  #
#=============#
### Created: ${date}
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: Getopt::Long${other_requirements}
#=============#
#  // HEADER  #
#=============#


#========================#
#  PRAGMAS & MODULES //  #
#========================#
use strict;
use warnings;
use Getopt::Long;
${other_pragmas_modules}
#========================#
#  // PRAGMAS & MODULES  #
#========================#


#======================#
#  USAGE & VERSION //  #
#======================#
sub usage {
### Returns usage information for current script in a string
<<USAGE;
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
sub version {
### Returns version information for current script in a string
<<VERSION;
$$0 version 1.0 (${date})
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on ${date}.
VERSION
}
#======================#
#  // USAGE & VERSION  #
#======================#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES <---#
my $$usage = 0;
my $$verbose = 0;
my $$version = 0;
${other_options_variables}

#---> PARSE / ASSIGN OPTIONS <---#
my $$options_result = GetOptions (
	${other_options}
	'usage|help' => \$$usage,
	'version' => \$$version,
	'verbose' => \$$verbose
);

#---> VERIFY OPTIONS <---#
# Print usage information and exit if option parsing was not successful
die &usage unless $$options_result;

# Print usage information and exit if --usage or --help were specified
die &usage if $$usage;

# Print version information and exit if --version was specified
die &version if $$version;

# Print usage information and exit if required arguments are not specified
die &usage if ! ${required_argument}; 
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#
#  MAIN //  #
#===========#
#---> STATUS MESSAGE <---#
print STDERR "Starting '$$0'...\n" if $$verbose;

#---> MAIN VARIABLES <---#
${main_variables}

	#---> BODY <---#

	#---> ${main_code_block_comment} <---#
	${main_code_block}

#---> STATUS MESSAGE <---#
print STDERR "Done.\n" if $$verbose;

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#


#==================#
#  SUBROUTINES //  #
#==================#
sub ${sub_name} {
### Function: ${sub_function_description}
### Accepts: ${sub_argument_description}
### Returns: ${sub_return_value_description}
### Dependencies: ${sub_dependencies}
### Type: ${sub_type_generic_specific}

	#---> PASS ARGUMENTS ---#
	${passed_variables}

	#---> STATUS MESSAGE <---#
	print STDERR "${sub_start_message}" . "\n" if $$verbose;

	#---> SUBROUTINE VARIABLES <---#
	${local_variables}

	#---> BODY <---#

		#---> ${sub_code_bloc} <---#
		${sub_code_block}

	#---> STATUS MESSAGE <---#
	print STDERR "${sub_end_message}" . "\n" if $$verbose;

	#---> RETURN VALUE <---#
	return ${return_value};

}
#==================#
#  // SUBROUTINES  #
#==================#