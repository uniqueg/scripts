#!/bin/sh

## Alexander Kanitz
## 14-MAR-2014

## Creates symbolic links of all files within the specified folder(s) (no recursion!) in $HOME/bin
## If no folder is specified, $PWD is used

# If no command line arguments are supplied, work on $PWD
if [ $# == 0 ]; then

        # Create symbolic links for each file
        for file in `find $PWD -maxdepth 1 -type f`; do
                base=`basename $file`
                ln -s $file $HOME/bin/$base
        done

# Else
else
        # ...traverse through all arguments passed on the command line
        for dir in $*; do

                # Skip if not a directory
                if [ -d $dir ]; then

                        # Get first character of directory
                        start=`echo $dir | cut -c 1,1`

                        # Add $PWD to directory if not an absolute path
                        # NOT ELEGANT FOR PATHES CONTAINING ".."
                        if [ ! "$start" == "/" ]; then
                                dir="$PWD/$dir"
                        fi

                        # Create symbolic links for each file
                        for file in `find $dir -maxdepth 1 -type f`; do
                                base=`basename $file`
                                ln -s $file $HOME/bin/$base
                        done

                else
                        echo -e "[WARNING] '$dir' is not a directory and will be ignored." >&2
                fi

        done

fi
~      