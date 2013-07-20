#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"
#include <omp.h>

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP: public PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType> {
	protected:
		size_t _numberOfThreads;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP(size_t maxRange, size_t blockSizeInElements = 16 * 1024, size_t numberOfThreads = 0, bool sendResultsToRoot = true,
				bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType>(maxRange, blockSizeInElements, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot), _numberOfThreads(
						numberOfThreads) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP() {
		}

		void removeMultiplesOfPrimesFromPreviousBlocksParallel(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t numberThreadsToUse = omp_get_max_threads();
			size_t numberOfThreads = this->template getNumberOfThreads();

			if (numberOfThreads != 0) {
				numberThreadsToUse = numberOfThreads;
			}

			size_t sievingMultiplesSize = sievingMultiples.size();

#			pragma omp parallel for \
				if (sievingMultiplesSize > 16) \
				default(shared) \
				schedule(guided) \
				num_threads(numberThreadsToUse)
			for (size_t sievingMultiplesIndex = 0; sievingMultiplesIndex < sievingMultiplesSize; ++sievingMultiplesIndex) {
				pair<size_t, size_t> primeCompositeInfo = sievingMultiples[sievingMultiplesIndex];
				size_t primeMultiple = primeCompositeInfo.first;
				size_t primeMultipleIncrement = primeCompositeInfo.second;

				for (; primeMultiple < blockEndNumber; primeMultiple += primeMultipleIncrement) {
					this->template setPrimesBitsetValueMPI(primeMultiple, true);
				}

				sievingMultiples[sievingMultiplesIndex].first = primeMultiple;
			}
		}


		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumberMPI(maxRangeSquareRoot);
			size_t blockBeginNumber = this->template getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t blockSizeInElements = this->template getBlockSizeInElements();

			size_t blockIndexBegin = this->template getBitsetPositionToNumberMPI(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + blockSizeInElements, maxIndexRangeSquareRoot + 1);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += blockSizeInElements;

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber >= maxRangeSquareRoot) {
					blockEndNumber = maxRangeSquareRoot + 1;
				}

//				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template removeMultiplesOfPrimesFromPreviousBlocksParallel(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}


		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
#			ifdef DEBUG_OUTPUT
			int _processID = this->template getProcessId();
			cout << "    --> Removing composites in process with rank " << _processID << " in [" << processBeginBlockNumber << ", " << (processEndBlockNumber - 1) << "]" << endl;
#			endif

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

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber > processEndBlockNumber) {
					blockEndNumber = processEndBlockNumber;
				}

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
#ifdef DEBUG_OUTPUT
				cout << "    > Probing for results..." << endl;
#endif
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);

				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
					int _numberProcesses = this->template getNumberProcesses();
					size_t processStartBlockNumber = this->template getProcessStartBlockNumber(processID, _numberProcesses, maxRange);
					if (processStartBlockNumber % 2 == 0) {
						++processStartBlockNumber;
					}

					size_t blockSize = this->template getProcessBitsetSize(processID, _numberProcesses, maxRange);
					size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(processStartBlockNumber);

					this->template receiveSievingDataMPI(primesBitset, positionToStoreResults, blockSize, status.MPI_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
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
				schedule(static) \
				reduction(+: primesFound) \
				num_threads(numberThreads)
			for (int threadBlockNumber = 0; threadBlockNumber < numberThreads; ++threadBlockNumber) {
				size_t possiblePrime;

				possiblePrime = threadBlockNumber * numberPrimesToCheckInBlock + startSieveNumber;
				if (!(wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t nextPossiblePrimeNumberEndBlock;
				if (threadBlockNumber == numberThreads - 1) {
					nextPossiblePrimeNumberEndBlock = maxRange + 1;
				} else {
					nextPossiblePrimeNumberEndBlock = (threadBlockNumber + 1) * numberPrimesToCheckInBlock + startSieveNumber;
				}

				while (possiblePrime < nextPossiblePrimeNumberEndBlock) {
					if (!this->template getPrimesBitsetValueMPI(possiblePrime)) {
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

