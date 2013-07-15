#include "PrimesSieveParallelMultiplesOptimizedOpenMPAndMPITimeAndCacheWithWheel.h"

template<>
void PrimesSieveParallelMultiplesOptimizedOpenMPAndMPITimeAndCacheWithWheel<PrimesFlagsContainerMPI, Modulo210WheelByte>::collectResultsFromProcessGroup(size_t maxRange) {
	cout << "\n    > Collecting results from other processes..." << endl;
	PrimesFlagsContainerMPI& primesBitset = getPrimesBitset();
	int numberProcesses = getNumberProcesses();

#pragma omp parallel for \
			default(shared) \
			firstprivate(numberProcesses, maxRange) \
			schedule(static)
	for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < numberProcesses; ++numberProcessesResultsCollected) {
		cout << "    > Probing for results..." << endl;
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD, &status);
		if (status.MPI_ERROR == MPI_SUCCESS) {
			int processID = status.MPI_SOURCE;
			size_t processStartBlockNumber = getProcessStartBlockNumber(processID, numberProcesses, maxRange);
			size_t processEndBlockNumber = getProcessEndBlockNumber(processID, numberProcesses, maxRange);

			if (processStartBlockNumber % 2 == 0) {
				++processStartBlockNumber;
			}

			if (processID == numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}
			size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;
			size_t positionToStoreResults = getBitsetPositionToNumberMPI(processStartBlockNumber);

			cout << "    --> Collecting results from process with rank " << processID << endl;
			MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD, &status);
			cout << "    --> Finished collecting results from process with rank " << processID << endl;
		} else {
			cout << "    --> MPI_Probe detected the following error code: " << status.MPI_ERROR << endl;
		}
	}
	cout << "    --> Finished collecting all results\n" << endl;
}
