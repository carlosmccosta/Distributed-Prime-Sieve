#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"

template<>
void PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>::sendResultsToRootProcess(size_t maxRange) {
	cout << "    > Sending results from process with rank " << _processID << " to root process..." << endl;
	PrimesFlagsContainerMPI& primesBitset = getPrimesBitset();
	size_t blockSize = getProcessBitsetSize(_processID, _numberProcesses, maxRange);

	MPI_Send(&primesBitset[0], blockSize, MPI_UNSIGNED_CHAR, 0, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD);
	cout << "    --> Finished sending results from process with rank " << _processID << " to root process" << endl;
}

template<>
void PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>::collectResultsFromProcessGroup(size_t maxRange) {
	cout << "\n    > Collecting results from other processes..." << endl;
	PrimesFlagsContainerMPI& primesBitset = getPrimesBitset();
	for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
		cout << "    > Probing for results..." << endl;
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK, MPI_COMM_WORLD, &status);
		if (status.MPI_ERROR == MPI_SUCCESS) {
			int processID = status.MPI_SOURCE;
			size_t processStartBlockNumber = getProcessStartBlockNumber(processID, _numberProcesses, maxRange);
			size_t blockSize = getProcessBitsetSize(_processID, _numberProcesses, maxRange);
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
