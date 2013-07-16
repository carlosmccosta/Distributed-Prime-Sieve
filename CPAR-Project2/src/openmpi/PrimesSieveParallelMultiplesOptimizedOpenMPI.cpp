#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"

template<>
MPI_Status PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>::receiveDataMPI(PrimesFlagsContainerMPI& primesBitset, size_t positionToStoreResults, size_t blockSize, int source,
		int tag) {
	MPI_Status status;
	cout << "    --> Collecting " << blockSize << " bytes of results from process with rank " << source << ", in process with rank " << _processID << endl;
	MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, source, tag, MPI_COMM_WORLD, &status);

	int numberBytesReceived;
	MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &numberBytesReceived);
	cout << "    --> Finished collecting " << numberBytesReceived << " bytes of results from process with rank " << source << ", in process with rank " << _processID << endl;
	return status;
}

template<>
void PrimesSieveParallelMultiplesOptimizedOpenMPI<PrimesFlagsContainerMPI, Modulo210WheelByte>::sendDataMPI(PrimesFlagsContainerMPI& primesBitset, size_t startPositionOfResults, size_t blockSize, int destination,
		int tag) {
	cout << "    --> Sending " << blockSize << " bytes of results from process with rank " << _processID << " to process with rank " << destination << endl;
	MPI_Send(&primesBitset[startPositionOfResults], blockSize, MPI_UNSIGNED_CHAR, destination, tag, MPI_COMM_WORLD);
	cout << "    --> Finished sending " << blockSize << " bytes of results from process with rank " << _processID << " to process with rank " << destination << endl;
}
