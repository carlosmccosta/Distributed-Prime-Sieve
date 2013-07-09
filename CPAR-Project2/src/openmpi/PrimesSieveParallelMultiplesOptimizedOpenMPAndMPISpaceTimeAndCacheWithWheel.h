#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel.h"

#include <omp.h>

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInElements = 16 * 1024) :
				PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements) {
		}
		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel() {
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = this->template getBlockSizeInElements();
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + processBeginBlockNumberIndex;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				if (blockIndexEnd > processEndBlockNumberIndex) {
					blockIndexEnd = processEndBlockNumberIndex;
				}

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockNumber == 0 && this->template getProcessId() == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
			cout << "\n    > Collecting results from other processes..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			MPI_Status status;
			int numProcesses = this->template getNumberProcesses();
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < numProcesses; ++numberProcessesResultsCollected) {
				cout << "    > Probing for results..." << endl;
				MPI_Probe(MPI_ANY_SOURCE, 1337, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, numProcesses, maxRange);
					size_t processEndBlockNumber = this->template getProcessEndBlockNumber(processID, numProcesses, maxRange);

					if (processStartBlockNumber % 2 == 0) {
						++processStartBlockNumber;
					}

					if (processID == numProcesses - 1) {
						processEndBlockNumber = maxRange + 1;
					}
					size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					cout << "    --> Collecting results from process with rank " << processID << endl;
					MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 1337, MPI_COMM_WORLD, &status);
					cout << "    --> Finished collecting results from process with rank " << processID << endl;
				} else {
					cout << "    --> MPI_Probe detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
			cout << "    --> Finished collecting all results\n" << endl;
		}

		virtual void sendResultsToRootProcess(size_t maxRange) {
			int processID = this->template getProcessId();
			int numberProcesses = this->template getNumberProcesses();
			cout << "    > Sending results from process with rank " << processID << " to root process..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, numberProcesses, maxRange);
			size_t processEndBlockNumber = this->template getProcessEndBlockNumber(processID, numberProcesses, maxRange);

			if (processStartBlockNumber % 2 == 0) {
				++processStartBlockNumber;
			}

			if (processID == numberProcesses - 1) {
				processEndBlockNumber = maxRange + 1;
			}
			size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;

			MPI_Send(&primesBitset[0], blockSize, MPI_UNSIGNED_CHAR, 0, 1337, MPI_COMM_WORLD);
			cout << "    --> Finished sending results from process with rank " << processID << " to root process" << endl;
		}

		virtual size_t getNumberPrimesFound() {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() >= 2)
				return primesValues.size();

			size_t primesFound = 4;
			size_t maxRange = this->template getMaxRange();
			int maxNumberThreads = omp_get_max_threads();
			size_t minNumberPrimesPerThread = 100;
			int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));
			size_t numberPrimesToCheckInBlock = maxRange / numberThreads;

			WheelType& wheelSieve = this->template getWheelSieve();

#pragma omp parallel for \
			default(shared) \
			firstprivate(maxRange, numberThreads, numberPrimesToCheckInBlock) \
			schedule(guided) \
			reduction(+: primesFound) \
			num_threads(numberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
				size_t possiblePrime;

				possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + 11;
				if (!(wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t nextPossiblePrimeNumberEndBlock = min((threadBlockNumber + 1) * numberPrimesToCheckInBlock + 11, maxRange + 1);

				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (this->template getPrimesBitsetValueMPI(possiblePrime)) {
						++primesFound;
					}
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}
			}

			return primesFound;
		}
};

