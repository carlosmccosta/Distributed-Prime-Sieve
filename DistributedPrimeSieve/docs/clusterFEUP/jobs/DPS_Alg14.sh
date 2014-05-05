#!/bin/bash
#PBS -j oe
#PBS -m abe
#PBS -M carlos.costa@fe.up.pt


echo "############################################################################################"
echo "# Modules"
echo "############################################################################################"
module load gcc/4.7.2
module load openmpi-x86_64
module list
echo -e "\n\n"


echo "############################################################################################"
echo "# PBS_NODEFILE"
echo "############################################################################################"
cat $PBS_NODEFILE
echo -e "\n\n"



echo "############################################################################################"
echo "# Program output"
echo "############################################################################################"
echo -e "\n\n"

mpiexec --prefix $MPI_HOME ~/DistributedPrimeSieve/Release/PrimeSieve --algorithm 14 --countPrimesInNode N --cacheBlockSize 16384 --segmentSizeInBlocks 65536 --numberThreads 16 --maxRangeInBits ${RANGE_IN_BITS}
