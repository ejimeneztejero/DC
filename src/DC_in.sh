#!/bin/bash

# This script is given to facilitate:
# - conversion of SGY/SEGY files into SU
# - and split it into smaller SU files of less than 1.9GB
# Perhaps this script has to be modified by the user requirements 

input_folder=$1	# path input data file
SGY_SU=$2 # file type
input=$3 # file name
NumRec=$4 #num channels
nt=$5 #time steps
ed_data=$6 #endianness data
ed_machine=$7 #endianness machine

TMPFILE="__check_header.su"

echo "Path input $key"  # path input data file
echo "File type $SGY_SU" # file type
echo "File name $input" # file name
echo "Num channels $NumRec" #num channels
echo "Time steps $nt" #time steps
if [ "$SGY_SU" = 1 ]; then
	if [ "$ed_data" = 0 ]; then
		echo "Endianness of data $ed_data" #endianness
	fi
	if [ "$ed_data" = 1 ]; then
		echo "Endianness of data $ed_data" #endianness
	fi
fi
echo "Endianness of machine $ed_machine" #endianness

cd $input_folder

echo "***** Finding if file is SGY or SU format *****"

echo "SGY_SU = $SGY_SU"

if surange < "$input" >/dev/null 2>&1; then
	echo "File seems SU format"
	SGY_SU=2
else
	echo "File is not SU format"
	SGY_SU=1
fi

echo "New SGY_SU = $SGY_SU"

echo "************************************************"

name=$(basename "$input" | cut -d. -f1)

if [ "$SGY_SU" = 1 ]; then

	echo "Convert .sgy file to .su file"

	su_name=${name}.su
	segyread tape="$input" conv=0 endian=0 segyclean > "$su_name"
	echo "segyread tape="$input" conv=0 endian=0 segyclean > "$su_name" "

fi	# sgy_su=1

if [ "$SGY_SU" = 2 ]; then
	su_name=${input}
fi

size=$(wc -c < "$su_name")
echo "size file: $size" ;

bytes1shot=$((NumRec*(nt+60)*4))
echo "bytes 1 shot: $bytes1shot";

maxbytes=1900000000

declare -i shots1file
shots1file=$((maxbytes/bytes1shot))
echo "shots in each part: $shots1file";

size_part=$((shots1file*bytes1shot))
echo "size (bytes): $size_part";

declare -i num_parts
quotient=$((size/size_part))
remainder=$((size%size_part))

if [ "$remainder" = 0 ]; then

	num_parts=$quotient

else

	num_parts=$((1+quotient))

fi

#echo "number of parts: $num_parts";

if [ $num_parts -le 10 ]
then
        num=1
fi

if [ $num_parts -gt 10 ] && [ $num_parts -le 100 ]
then
        num=2
fi

if [ $num_parts -gt 100 ] && [ $num_parts -le 1000 ]
then
        num=3
fi


if [ $num_parts = 1 ]
then

echo
echo "No need to split SU file, it has already correct size"
echo

name_part=$su_name

fi

if [ $num_parts -gt 1 ]
then

echo
echo "Split SU file in SU files of less than 1.9 GB"
echo

name_part=${name}_

####################################################################

# Perhaps the following line has to be modified by the user:
echo "split -b $size_part -d -a $num $su_name $name_part"
split -b $size_part -d -a $num $su_name $name_part

####################################################################

fi

rm DC_in.dat
echo "$name_part" > DC_in.dat
echo "$num_parts" >> DC_in.dat
#echo "$ed_data" >> DC_in.dat

rm binary
rm header

####################################################################
