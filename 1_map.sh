#!/bin/bash
#Made by Mykhaylo Usyk MSci. 2016 mu408@nyu.edu

#Reference Directories and Files
main_dir=$(pwd)
mast="${1}"

#Set fancy fonts for the help message
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

#Help
function HELP {
  echo "${BOLD}Help documentation for the Map generator${NORM}"
  echo "The following options must be specified"
  echo "${REV}-m${NORM}   The location of the ${BOLD}pcr text file${NORM} that is a tab delimited text with Unix line breaks"
  echo "The following options can be passed to the script, but are not mandatory:"
  echo "${REV}-h${NORM}   Displays this extremely helpful message and aborts script"
  echo "Example: sh map.sh${BOLD} -m /home/Desktop/pcr_map.txt"
  exit 1
}

#Check the number of arguments. If none are passed, print message and exit
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  echo "Did not pass any arguments, you probably need help\n"
  HELP
fi

#Parse the inputs
while getopts :m:h FLAG; do
  case $FLAG in
    m)  #set option "m"
      OPT_m=$OPTARG
      ;;
    h)  #set option "h"
      OPT_p=$OPTARG
      HELP
      ;;
    \?) #unrecognized option - show help
      echo "Option -${BOLD}$OPTARG${NORM} not allowed."
      exit 1
      ;;
  esac
done



if [[ -z "$OPT_m" ]]; then
	echo "No map file given, aborting script"
	exit 1
fi


#get the plate barcodes for the genotype

cat ${OPT_m} | cut -f 4,5 | sed '1d' | sort | uniq | tr '\t' ' ' > 1_1_1_plate_barcodes.txt

# Make general barcodes
echo "Distance 1" > general_barcode_p1.txt
echo "Format 5" >> general_barcode_p1.txt

cat 1_1_1_plate_barcodes.txt >> general_barcode_p1.txt

echo "Distance 1" > general_barcode_p2.txt
echo "Format N 5" >> general_barcode_p2.txt

cat 1_1_1_plate_barcodes.txt >> general_barcode_p2.txt

echo "Distance 1" > 1_1_3_forward.txt
echo "Format 5" >> 1_1_3_forward.txt

cat ${OPT_m} | cut -f 2,3 | sed '1d' | sort | uniq | tr '\t' ' ' >> 1_1_3_forward.txt