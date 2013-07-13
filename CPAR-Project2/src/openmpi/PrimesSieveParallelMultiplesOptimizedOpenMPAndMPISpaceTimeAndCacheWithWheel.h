#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel.h"

#include <omp.h>

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType> {
	protected:
		size_t _numberOfThreads;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel<FlagsContainer, WheelType>(maxRange, blockSizeInElements, sendResultsToRoot, sendPrimesCountToRoot), _numberOfThreads(numberOfThreads) {
		}
		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPAndMPISpaceTimeAndCacheWithWheel() {
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = this->template getBlockSizeInElements();
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);
			size_t numberThreadsToUse = omp_get_max_threads();
			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					numberThreadsToUse = numberBlocks;
				} else {
					numberThreadsToUse = _numberOfThreads;
				}
			}

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

#pragma omp parallel for \
			default(shared) \
			firstprivate(processBeginBlockNumberIndex, processEndBlockNumberIndex, blockSizeInElements, sievingMultiples, priviousBlockNumber) \
			schedule(guided, 64) \
			num_threads(numberThreadsToUse)
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
			int numberProcesses = this->template getNumberProcesses();

#pragma omp parallel for \
			default(shared) \
			firstprivate(numberProcesses, maxRange) \
			schedule(static)
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < numberProcesses; ++numberProcessesResultsCollected) {
				cout << "    > Probing for results..." << endl;
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, numberProcesses, maxRange);
					size_t processEndBlockNumber = this->template getProcessEndBlockNumber(processID, numberProcesses, maxRange);

					if (processStartBlockNumber % 2 == 0) {
						++processStartBlockNumber;
					}

					if (processID == numberProcesses - 1) {
						processEndBlockNumber = maxRange + 1;
					}
					size_t blockSize = ((processEndBlockNumber - processStartBlockNumber) >> 1) + 1;
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					cout << "    --> Collecting results from process with rank " << processID << endl;
					MPI_Recv(&(primesBitset[positionToStoreResults]), blockSize, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
					cout << "    --> Finished collecting results from process with rank " << processID << endl;
				} else {
					cout << "    --> MPI_Probe detected the following error code: " << status.MPI_ERROR << endl;
				}
			}
			cout << "    --> Finished collecting all results\n" << endl;
		}


		virtual size_t getNumberPrimesFound() {
			size_t primesFound = this->template getPrimesCount();
			// avoid recomputation
			if (primesFound != 0) {
				return primesFound;
			}

			vector<size_t>& primesValues = this->template getPrimesValues();
			if (primesValues.size() >= 2) {
				return primesValues.size();
			}

			size_t startSieveNumber = this->template getStartSieveNumber();
			if (startSieveNumber == this->template getBlockBeginNumber()) {
				primesFound = this->template getWheelSieve().getNumberPrimesSievedByTheWheel();
			} else {
				primesFound = 0;
			}

			size_t maxRange = this->template getMaxRange();
			int maxNumberThreads = omp_get_max_threads();
			size_t minNumberPrimesPerThread = 100;
			int numberThreads = min((size_t) maxNumberThreads, (size_t) ceil((double) maxRange / (double) minNumberPrimesPerThread));
			size_t numberPrimesToCheckInBlock = (maxRange - startSieveNumber) / numberThreads;

			WheelType& wheelSieve = this->template getWheelSieve();

#pragma omp parallel for \
			default(shared) \
			firstprivate(maxRange, numberThreads, numberPrimesToCheckInBlock, startSieveNumber) \
			schedule(guided) \
			reduction(+: primesFound) \
			num_threads(numberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
				size_t possiblePrime;

				possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + startSieveNumber;
				if (!(wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t nextPossiblePrimeNumberEndBlock = min((threadBlockNumber + 1) * numberPrimesToCheckInBlock + startSieveNumber, maxRange + 1);

				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (this->template getPrimesBitsetValueMPI(possiblePrime)) {
						++primesFound;
					}
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}
			}

			this->template setPrimesCount(primesFound);
			return primesFound;
		}

		size_t getNumberOfThreads() const {
			return _numberOfThreads;
		}

		void setNumberOfThreads(size_t numberOfThreads) {
			_numberOfThreads = numberOfThreads;
		}
};

