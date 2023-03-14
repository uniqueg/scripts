#!/bin/bash

# execute command of which the exit status should be checked

exit_status=$?

if [ $exit_status -eq 0 ]
then
	echo "Exit status: $exit_status"
fi

if [ $exit_status -eq 1 ]
then
	echo "Exit status: $exit_status"
fi