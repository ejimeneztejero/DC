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

sgy_DC="DC0.sgy"
su_DC="DC0.su"

echo
echo 'DC0:'
echo "First step: join all files in one"

echo
read -p "Name of files without the number: " DC_part
cat ${DC_part}* > ${su_DC}

echo "Second step: convert .su file to .sgy file"

segyhdrs < $su_DC
segywrite < $su_DC tape=${sgy_DC}
echo ''

fi

####################################################################

if [[ $answer1 = yes ]]; then 

su_DC="DC1.su"
sgy_DC="DC1.sgy"

echo
echo 'DC1'

echo "First step: join all files in one"

echo
read -p "Name of files without the number: " DC_part
cat ${DC_part}* > ${su_DC}

echo "Second step: convert .su file to .sgy file"
segyhdrs < $su_DC
segywrite < $su_DC tape=${sgy_DC}
echo ''

fi

####################################################################

if [[ $answer2 = yes ]]; then 

su_DC="DC2.su"
sgy_DC="DC2.sgy"

echo
echo 'DC2'

echo "First step: join all files in one"

echo
read -p "Name of files without the number: " DC_part
cat ${DC_part}* > ${su_DC}

echo "Second step: convert .su file to .sgy file"
segyhdrs < $su_DC
segywrite < $su_DC tape=${sgy_DC}
echo ''

fi

####################################################################
