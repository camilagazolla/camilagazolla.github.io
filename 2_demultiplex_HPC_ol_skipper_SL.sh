#!/bin/bash
#Made by Mykhaylo Usyk MSci. 2022 mu408@nyu.edu

# module load novocraft/Oct20
PATH=/gs/gsfs0/users/burk-lab/Biosoftware/:$PATH
# module load prinseq-lite/0.20.4
PATH=/gs/gsfs0/users/burk-lab/BioSoftware_CZ/prinseq-lite-0.20.4/:$PATH

#Important directories/files
main_dir=$(pwd)

#Novobarcode path
#novo="/Users/LisaK/Desktop/ana/ana_test/BioSoftware/novocraft/novobarcode"
#prinseq="/Users/LisaK/Desktop/ana/ana_test/BioSoftware/prinseq-lite-0.20.4/prinseq-lite.pl"

#Set fancy fonts for the help message
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

#Help
function HELP {
  echo "${BOLD}Help documentation for the Demultiplexer${NORM}"
  echo "Please specify the forward and reverse read locations"
  echo "Make sure to have a fastq directory with the forward and reverse reads in the same location as this script"
  echo "${REV}-f${NORM}   The location of the ${BOLD}forward reads${NORM}"
  echo "${REV}-r${NORM}   The location of the ${BOLD}reverse reads${NORM}"
  echo "${REV}-m${NORM}   The location of the ${BOLD}original mapping file${NORM}"
  echo "${REV}-s${NORM}   Passing this option skips the prinseq-lite step"
  echo "${REV}-h${NORM}   Displays this extremely helpful message and aborts script"
  echo "Example: sh general_demultiplex.sh${BOLD} -f fastq/forward_R1.fastq -r fastq/reverse_R1.fastq${NORM}"
  exit 1
}

#Check the number of arguments. If none are passed, print message and exit
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  echo "Did not pass any arguments, you probably need help\n"
  HELP
fi

# define the run type of this script
runtype="Everything"

#Parse the inputs
while getopts :f:r:m:hs FLAG; do
  case $FLAG in
    f)  #set option "f"
      OPT_f=$OPTARG
      ;;
    r)  #set option "r"
      OPT_r=$OPTARG
      ;;
    m)  #set option "m"
      OPT_m=$OPTARG
      ;;
    h)  #set option "h"
      OPT_p=$OPTARG
      HELP
      ;;
    s)  # prinseq skippy thing 
      runtype="PrinseqSkip"
      ;;
    \?) #unrecognized option - show help
      echo "Option -${BOLD}$OPTARG${NORM} not allowed."
      exit 1
      ;;
  esac
done

if [[ -z "$OPT_f" ]]; then
	echo "No forward read specified, aborting script"
	exit 1
fi

if [[ -z "$OPT_r" ]]; then
	echo "No reverse read specified, aborting script"
	exit 1
fi

if [[ -z "$OPT_m" ]]; then
	echo "No mapping file, aborting script"
	exit 1
fi


#Library variable
lib="libA"
mkdir ${main_dir}/fastq
cd ${main_dir}/fastq

if [ $runtype = "Everything" ]; then

echo "Full Pipe Selected, running prinseq-lite"

#Trimming the 3 bp buffers
prinseq-lite \
-fastq ${OPT_f} \
-fastq2 ${OPT_r} \
-out_format 3 \
-out_bad null \
-trim_left 3 \
-no_qual_header \
-out_good p29_${lib} \
-trim_qual_right 25 \
-trim_qual_left 25 \
-min_len 50 \
-verbose

# just doing a simple else without any kind of check, because I won't be adding any other runtype options (hopefully)
else 

echo "Skipping prinseq-lite step"

fi


cd ..

# create two folders for general reverse barcode demultiplexing #
mkdir -p ${main_dir}/novobarcode/p29_General_Barcode_Golay_P1/
mkdir -p ${main_dir}/novobarcode/p29_General_Barcode_Golay_P2/
mkdir -p ${main_dir}/novobarcode/log/


# demultiplex general (reverse) barcoded reads #
novobarcode_HPC \
-b ${main_dir}/general_barcode_p1.txt \
-d ${main_dir}/novobarcode/p29_General_Barcode_Golay_P1 \
-f ${main_dir}/fastq/p29_${lib}_1.fastq ${main_dir}/fastq/p29_${lib}_2.fastq \
--NC_OFF \
> ${main_dir}/novobarcode/log/p29_General_Barcode_Golay_P1_log.txt

novobarcode_HPC \
-b ${main_dir}/general_barcode_p2.txt \
-d ${main_dir}/novobarcode/p29_General_Barcode_Golay_P2 \
-f ${main_dir}/fastq/p29_${lib}_1.fastq ${main_dir}/fastq/p29_${lib}_2.fastq \
--NC_OFF \
> ${main_dir}/novobarcode/log/p29_General_Barcode_Golay_P2_log.txt

echo "Finished General Demultiplex"

cd ${main_dir}

#Library
lib="libA"

ReverseListPAP="${main_dir}/1_1_1_plate_barcodes.txt"

UniqueBCPAP="${main_dir}/1_1_3_forward.txt"
  

line_no=$(wc -l < $ReverseListPAP)
  
for i in $(seq 1 ${line_no}); do \
R_bc_list=$(sed "${i}q;d" $ReverseListPAP); \
Reverse_List_Name=($(echo "${R_bc_list}" | cut -d " " -f 1));\
Reverse_List_Seq=($(echo "${R_bc_list}" | cut -d " " -f 2)); \
echo "Processing ${Reverse_List_Name} demultiplexing ..."; \
mkdir -p ${main_dir}/novobarcode/${Reverse_List_Name}/; \
mkdir -p ${main_dir}/novobarcode/concatinated_reads_${lib}/; \
wait; \
cat ${main_dir}/novobarcode/p29_General_Barcode_Golay_P1/${Reverse_List_Seq}/p29_${lib}_1.fastq \
${main_dir}/novobarcode/p29_General_Barcode_Golay_P2/${Reverse_List_Seq}/p29_${lib}_2.fastq > \
${main_dir}/novobarcode/concatinated_reads_${lib}/p29_${lib}_${Reverse_List_Name}_P2.fastq; \
wait; \
cat ${main_dir}/novobarcode/p29_General_Barcode_Golay_P1/${Reverse_List_Seq}/p29_${lib}_2.fastq \
${main_dir}/novobarcode/p29_General_Barcode_Golay_P2/${Reverse_List_Seq}/p29_${lib}_1.fastq > \
${main_dir}/novobarcode/concatinated_reads_${lib}/p29_${lib}_${Reverse_List_Name}_P1.fastq; \
wait; \
novobarcode_HPC \
-b ${UniqueBCPAP} \
-d ${main_dir}/novobarcode/${Reverse_List_Name}/ \
-f ${main_dir}/novobarcode/concatinated_reads_${lib}/p29_${lib}_${Reverse_List_Name}_P1.fastq \
${main_dir}/novobarcode/concatinated_reads_${lib}/p29_${lib}_${Reverse_List_Name}_P2.fastq > \
${main_dir}/novobarcode/log/p29_${lib}_loop_${Reverse_List_Name}_log.txt; \
done

echo "Cleaning Up"

cd ${main_dir}/novobarcode
mkdir demultiplex_loop_1
mv Golay* demultiplex_loop_1/

echo "Done Forward Demultiplex"

cd ${main_dir}

sed '1d' ${OPT_m} > to_get_samples

mkdir samples

while read p
do 
sampleid=$(echo "${p}" | cut -f 1)
rev_bar=$(echo "${p}" | cut -f 4)
for_seq=$(echo "${p}" | cut -f 3)
mkdir samples/${sampleid}
cp novobarcode/demultiplex_loop_1/${rev_bar}/${for_seq}/* samples/${sampleid}
done < to_get_samples

cd samples

for f in */; do sampid=$(echo "${f}" | sed 's|/||1')
mv ${f}*P1.fastq ${f}${sampid}_R1.fastq
mv ${f}*P2.fastq ${f}${sampid}_R2.fastq
done

cd ${main_dir}

echo "Finished moving samples"





