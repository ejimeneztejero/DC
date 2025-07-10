#!/bin/bash

output_folder=$1	# path file
file_part=$2 		# file name
num=$3			# DC file, 0, 1 or 2

cd $output_folder

####################################################################

sgy_DC="DC$num.sgy"
su_DC="DC$num.su"

echo
echo "First step: join all files ($file_part) into one ($su_DC)"

cat ${file_part}* > ${su_DC}

echo "Second step: convert $su_DC file to $sgy_DC file"

segyhdrs < $su_DC
segywrite < $su_DC tape=${sgy_DC}
echo

mv $sgy_DC ..
mv $su_DC ..

rm binary
rm header

####################################################################
