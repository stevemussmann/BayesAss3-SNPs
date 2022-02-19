#!/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
	# Display Help
	echo "This script will report the number of loci in"
	echo "an immanc-formatted file."
	echo ""
	echo "The following are the options for this script."
	echo ""
	echo "Syntax: countLociImmanc.sh [-f|h]"
	echo "options:"
	echo "f		Use this to input an immanc file."
	echo "h		Prints this help menu."
	echo ""
}

# parse the command line options.
while getopts "f:h" option; do
	case $option in
		h) # display Help
			Help
			exit;;
		f) # set file name
			FILE=$OPTARG
	esac
done

# check if there were any command line options provided and exit if none detected
if [ $OPTIND -eq 1 ];
then
		echo ""
		echo "No options were passed."
		echo "Use -h to display the help menu."
		echo ""
		exit
fi

# count loci in immanc file

if [ -f $FILE ]
then
	echo ""
	echo $FILE
	awk '{print $3}' $FILE | sort | uniq | wc -l
	echo ""
else
	echo ""
	echo "$FILE does not exist."
	echo ""
fi


exit
