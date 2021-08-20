#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -e make_error
#$ -o make_output
#$ -q new_ifb.q

module purge
module load intel-2017U2 openmpi-2.1.1_ifb-intel2017U2

export PATH=/opt/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi/lib:/usr/lib64/libblas.so.3:$LD_LIBRARY_PATH

make


