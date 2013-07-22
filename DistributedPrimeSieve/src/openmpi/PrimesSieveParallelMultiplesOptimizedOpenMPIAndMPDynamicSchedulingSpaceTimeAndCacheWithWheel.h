#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, size_t numberOfThreads = 0,
				bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true, bool sendPrimesCountToRoot = true, size_t dynamicSchedulingSegmentSizeInElements = 1048576,
				size_t dynamicSchedulingNumberSegments = 0, string outputResultsFilename = "", bool outputOnlyLastSegment = false) :
				PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, numberOfThreads, sendResultsToRoot,
						countNumberOfPrimesOnNode, sendPrimesCountToRoot, dynamicSchedulingSegmentSizeInElements, dynamicSchedulingNumberSegments, outputResultsFilename, outputOnlyLastSegment) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicSchedulingSpaceTimeAndCacheWithWheel() {
		}

		virtual inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			if (blockSize % 2 == 0) {
				return (blockSize >> 1);
			} else {
				return ((blockSize >> 1) + 1);
			}
		}

		inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
			return (position << 1) + this->template getStartSieveNumber();
		}

		// ------------------------------ <code from PrimesSieveParallelMultiplesOptimizedOpenMPI to avoid using virtuals in memory access functions (hotspots)> -------------------------------------------

		inline bool getPrimesBitsetValueMPI(size_t number) {
			return this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)];
		}

		inline void setPrimesBitsetValueMPI(size_t number, bool newValue) {
			this->template getPrimesBitset()[this->template getBitsetPositionToNumberMPI(number)] = newValue;
		}

		void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t sievingMultiplesSize = sievingMultiples.size();
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

		void calculatePrimesInBlock(size_t blockBeginNumber, size_t blockEndNumber, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxPrimeNumberSearch = blockEndNumber;
			if (maxPrimeNumberSearch >= maxRangeSquareRoot) {
				maxPrimeNumberSearch = maxRangeSquareRoot + 1;
			}

			size_t primeNumber = blockBeginNumber;
			WheelType _wheelSieve = this->template getWheelSieve();
			if (_wheelSieve.getBitsetPositionToNumberWithCheck(primeNumber) == std::numeric_limits<std::size_t>::max()) {
				primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber);
			}

			vector<size_t>& _sievingPrimes = this->template getSievingPrimes();

			for (; primeNumber < maxPrimeNumberSearch; primeNumber = _wheelSieve.getNextPossiblePrime(primeNumber)) {
				// for each number not marked as composite (prime number)
				if (!this->template getPrimesBitsetValueMPI(primeNumber)) {
					_sievingPrimes.push_back(primeNumber);

					//use it to calculate his composites
					size_t primeMultipleIncrement = primeNumber << 1;
					size_t compositeNumber = primeNumber * primeNumber;
					for (; compositeNumber < blockEndNumber; compositeNumber += primeMultipleIncrement) {
						this->template setPrimesBitsetValueMPI(compositeNumber, true);
					}
					sievingMultiples.push_back(pair<size_t, size_t>(compositeNumber, primeMultipleIncrement));
				}
			}
		}

		virtual void savePrimes(ostream& outputStream) {
			vector<size_t>& primesValues = this->template getPrimesValues();

			if (primesValues.size() <= 2) {
				size_t possiblePrime = this->template getStartSieveNumber();

				if (possiblePrime == this->template getWheelSieve().getFirstPrimeToSieve()) {
					outputStream << 2 << endl;
					outputStream << 3 << endl;
					outputStream << 5 << endl;
					outputStream << 7 << endl;
				}

				WheelType _wheelSieve = this->template getWheelSieve();

				if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
					possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
				}

				size_t maxRange = this->template getMaxRange();
				for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
					if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
						outputStream << possiblePrime << endl;
					}
				}
			} else {
				size_t iSize = primesValues.size();
				for (size_t i = 0; i < iSize; ++i) {
					outputStream << primesValues[i] << endl;
				}
			}
		}

		virtual vector<size_t>& extractPrimesFromBitset() {
			vector<size_t>& primesValues = this->template getPrimesValues();

			primesValues.clear();
			primesValues.push_back(2);
			primesValues.push_back(3);
			primesValues.push_back(5);
			primesValues.push_back(7);

			size_t possiblePrime = this->template getStartSieveNumber();
			size_t maxRange = this->template getMaxRange();
			WheelType _wheelSieve = this->template getWheelSieve();

			for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
		}
		// ------------------------------ </code from PrimesSieveParallelMultiplesOptimizedOpenMPI to avoid using virtuals in memory access functions (hotspots)> ------------------------------------------

		// ---------------------------- <code from PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP to avoid using virtuals in memory access functions (hotspots)> ----------------------------------------
		virtual void removeMultiplesOfPrimesFromPreviousBlocksParallel(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t numberThreadsToUse = omp_get_max_threads();
			size_t numberOfThreads = this->template getNumberOfThreads();

			if (numberOfThreads != 0) {
				numberThreadsToUse = numberOfThreads;
			}

			size_t sievingMultiplesSize = sievingMultiples.size();

#			pragma omp parallel for \
						if (sievingMultiplesSize > 32) \
						default(shared) \
						schedule(guided, 8) \
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
			size_t _numberOfThreads = this->template getNumberOfThreads();
			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					numberThreadsToUse = numberBlocks;
				} else {
					numberThreadsToUse = _numberOfThreads;
				}
			}

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

#			pragma omp parallel for \
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

#			pragma omp parallel for \
						default(shared) \
						firstprivate(numberProcesses, maxRange) \
						schedule(static)
			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < numberProcesses; ++numberProcessesResultsCollected) {
#				ifdef DEBUG_OUTPUT
				cout << "    > Probing for results..." << endl;
#				endif
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

#			pragma omp parallel for \
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
		// ---------------------------- </code from PrimesSieveParallelMultiplesOptimizedOpenMPIAndMP to avoid using virtuals in memory access functions (hotspots)> ---------------------------------------

		// ------------------------<code from PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling to avoid using virtuals in memory access functions (hotspots)> ----------------------------
		size_t receiveSievedBlockFromNode() {
			// received block range
			unsigned long long segmentRange[2];
			MPI_Status status;

#			ifdef DEBUG_OUTPUT
			cout << "    --> Ready to receive segment results from sieving node..." << endl;
#			endif

			MPI_Recv(&segmentRange[0], 2, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_BLOCK_RANGE, MPI_COMM_WORLD, &status);

			if (status.MPI_ERROR == MPI_SUCCESS) {
				size_t blockSize = this->template getNumberBitsToStoreBlock(segmentRange[1] - segmentRange[0]);
				size_t positionToStoreResults = this->template getBitsetPositionToNumberMPI(segmentRange[0]);

				// receive sieving data
				FlagsContainer& primesBitset = this->template getPrimesBitset();
				this->template receiveSievingDataMPI(primesBitset, positionToStoreResults, blockSize, status.MPI_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT);
				return blockSize;
			} else {
				cerr << "    !!!!! Detected error " << status.MPI_ERROR << " when receiving sieving segment results data from process with rank " << status.MPI_SOURCE << "!!!!!" << endl;
				return 0;
			}
		}
		// -----------------------</code from PrimesSieveParallelMultiplesOptimizedOpenMPIAndMPDynamicScheduling to avoid using virtuals in memory access functions (hotspots)> ----------------------------

};

