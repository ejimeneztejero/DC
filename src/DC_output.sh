#!/bin/bash

# Perhaps the following lines has to be modified by the user:

echo
read -p 'Press enter if you are already in output folder or write path (ex. /home/user/folder/) :' key
echo
read -p "DC0: join files and convert to .sgy [yes or no]: " answer0
echo
read -p "DC1: join files and convert to .sgy [yes or no]: " answer1
echo
read -p "DC2: join files and convert to .sgy [yes or no]: " answer2

if [ "$key" = '' ]; then
        output_folder=${PWD}
else
    	output_folder=$key
fi

cd $output_folder

####################################################################

if [[ $answer0 = yes ]]; then 

echo
echo 'DC0:'
echo "First step: join all files in one"

file_DC0=su_output_DC0
DC0_part=${file_DC0}_part_
cat ${DC0_part}* > ${file_DC0}

echo "Second step: convert .su file to .sgy file"

sgy_DC0=output_DC0
segyhdrs < $file_DC0
segywrite < $file_DC0 tape=${sgy_DC0}.sgy
echo ''

fi

####################################################################

if [[ $answer1 = yes ]]; then 

echo
echo 'DC1'

echo "First step: join all files in one"
file_DC1=su_output_DC1
DC1_part=${file_DC1}_part_
cat ${DC1_part}* > ${file_DC1}

echo "Second step: convert .su file to .sgy file"
sgy_DC1=output_DC1
segyhdrs < $file_DC1
segywrite < $file_DC1 tape=${sgy_DC1}.sgy
echo ''

fi

####################################################################

if [[ $answer2 = yes ]]; then 

echo
echo 'DC2'

echo "First step: join all files in one"
file_DC2=su_output_DC2
DC2_part=${file_DC2}_part_
cat ${DC2_part}* > ${file_DC2}

echo "Second step: convert .su file to .sgy file"
sgy_DC2=output_DC2
segyhdrs < $file_DC2
segywrite < $file_DC2 tape=${sgy_DC2}.sgy

fi

####################################################################
