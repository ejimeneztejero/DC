#!/bin/bash

echo
read -p 'Press enter if you are already in input folder or write path (ex. /home/user/folder/) :' key
read -p 'Are you reading SGY file (press 1) or SU file (press 2): ' SGY_SU
read -p 'Name of your file : ' input
read -p 'Number of channels: ' NumRec
read -p 'ns: ' nt
echo

if [ "$key" = '' ]; then
        input_folder=${PWD}
else
        input_folder=$key
fi

cd $input_folder


name=$(basename "$input" | cut -d. -f1)
su_name=${name}


if [ "$SGY_SU" = 1 ]; then


	echo "Convert .sgy file to .su file"
	echo

	# Perhaps the following line has to be modified by the user:

	segyread tape=$input verbose=0 conv=0 endian=0 segyclean > $su_name

fi 


name_part=${su_name}_part_

size=$(wc -c < "$su_name")
echo "size file: $size" ;

bytes1shot=$((NumRec*(nt+60)*4))
echo "bytes 1 shot: $bytes1shot";

maxbytes=1900000000

declare -i shots1file
shots1file=$((maxbytes/bytes1shot))
echo "shots in each part: $shots1file";

size_part=$((shots1file*bytes1shot))
echo "size (bytes) in each part: $size_part";

declare -i num_parts
quotient=$((size/size_part))
remainder=$((size%size_part))

if [ "$remainder" = 0 ]; then

	num_parts=$quotient

else

	num_parts=$((1+quotient))

fi

echo "number of parts: $num_parts";

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

fi

if [ $num_parts -gt 1 ]
then

echo
echo "Split SU file in SU files of less than 1.9 GB"
echo

####################################################################

# Perhaps the following line has to be modified by the user:
split -b $size_part -d -a $num $su_name $name_part

####################################################################

fi
