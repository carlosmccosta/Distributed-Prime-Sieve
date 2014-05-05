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


# Filter PBS_NODEFILE to remove repeated node hosts 
nodesFile=$PBS_O_WORKDIR/$PBS_JOBID.nodes
sort $PBS_NODEFILE | uniq > $nodesFile



echo "############################################################################################"
echo "# PBS_NODEFILE filtered"
echo "############################################################################################"
cat $nodesFile
echo -e "\n\n"


# assign the specified number of slots to each node host
hostsFileProcessed=$PBS_O_WORKDIR/$PBS_JOBID.hostfile
for jobNode in $(cat $nodesFile); do
	echo $jobNode slots=$NUMBER_SLOTS_PER_NODE >> $hostsFileProcessed
done



echo "############################################################################################"
echo "# Final hostfile"
echo "############################################################################################"
cat $hostsFileProcessed
echo -e "\n\n"



echo "############################################################################################"
echo "# Program output"
echo "############################################################################################"
echo -e "\n\n"

mpiexec --prefix $MPI_HOME --np $NUMBER_PROCESSES --hostfile $hostsFileProcessed ~/DistributedPrimeSieve/Release_2/PrimeSieve --algorithm 19 --countPrimesInNode N --sendPrimesCountToRoot N --sendResultsToRoot N --outputOnlyLastSegment N --cacheBlockSize 16384 --numberThreads 8 --dynamicSchedulingNumberSegments 4096 --maxRangeInBits ${RANGE_IN_BITS}

rm -f $nodesFile
rm -f $hostsFileProcessed
