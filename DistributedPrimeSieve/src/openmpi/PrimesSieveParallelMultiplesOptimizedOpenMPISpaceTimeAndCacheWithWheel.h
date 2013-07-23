#pragma once

#include "PrimesSieveParallelMultiplesOptimizedOpenMPI.h"

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel: public PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType> {
	public:
		PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel(size_t maxRange, size_t blockSizeInBytes = 16 * 1024, bool sendResultsToRoot = true, bool countNumberOfPrimesOnNode = true,
				bool sendPrimesCountToRoot = true) :
				PrimesSieveParallelMultiplesOptimizedOpenMPI<FlagsContainer, WheelType>(maxRange, blockSizeInBytes * 8, sendResultsToRoot, countNumberOfPrimesOnNode, sendPrimesCountToRoot) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMPISpaceTimeAndCacheWithWheel() {
		}

		virtual inline size_t getNumberBitsToStoreBlock(size_t blockSize) {
			if (blockSize % 2 == 0) {
				return (blockSize >> 1);
			} else {
				return ((blockSize >> 1) + 1);
			}
		}

		virtual inline size_t getBitsetPositionToNumberMPI(size_t number) {
			return (number - this->template getStartSieveNumber()) >> 1;
		}

		virtual inline size_t getNumberAssociatedWithBitsetPositionMPI(size_t position) {
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


		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumberMPI(maxRangeSquareRoot);
			size_t blockBeginNumber = this->template getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t _blockSizeInElements = this->template getBlockSizeInElements();
			size_t blockIndexBegin = this->template getBitsetPositionToNumberMPI(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + _blockSizeInElements, maxIndexRangeSquareRoot + 1);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

			size_t numberBlocks = (size_t) ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber >= maxRangeSquareRoot) {
					blockEndNumber = maxRangeSquareRoot + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		virtual void removeComposites(size_t processBeginBlockNumber, size_t processEndBlockNumber, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
#			ifdef DEBUG_OUTPUT
			cout << "    --> Removing composites in process with rank " << _processID << " in [" << processBeginBlockNumber << ", " << (processEndBlockNumber - 1) << "]" << endl;
#			endif

			size_t _blockSizeInElements = this->template getBlockSizeInElements();
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t processEndBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processEndBlockNumber);
			const size_t processBeginBlockNumberIndex = this->template getBitsetPositionToNumberMPI(processBeginBlockNumber);

			const size_t numberBlocks = (size_t) ceil((double) (processEndBlockNumberIndex - processBeginBlockNumberIndex) / (double) blockSizeInElements);

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockEndNumber = -1;

			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + processBeginBlockNumberIndex;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionMPI(blockIndexEnd);

				if (blockEndNumber > processEndBlockNumber) {
					blockEndNumber = processEndBlockNumber;
				}

				size_t _processID = this->template getProcessId();

				if (blockNumber == 0 && _processID == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockBeginNumber != priviousBlockEndNumber) {
					this->template computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockEndNumber = blockEndNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void collectResultsFromProcessGroup(size_t maxRange) {
			cout << "\n    > Collecting results from other processes..." << endl;
			FlagsContainer& primesBitset = this->template getPrimesBitset();
			int _numberProcesses = this->template getNumberProcesses();

			for (int numberProcessesResultsCollected = 1; numberProcessesResultsCollected < _numberProcesses; ++numberProcessesResultsCollected) {
#				ifdef DEBUG_OUTPUT
				cout << "    > Probing for results..." << endl;
#				endif
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MSG_NODE_COMPUTATION_RESULTS_SEGMENT, MPI_COMM_WORLD, &status);
				if (status.MPI_ERROR == MPI_SUCCESS) {
					int processID = status.MPI_SOURCE;
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

			size_t maxRange = this->template getMaxRange();
			WheelType _wheelSieve = this->template getWheelSieve();

			size_t possiblePrime = this->template getStartSieveNumber();
			if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
				possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
			}

			for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
					primesValues.push_back(possiblePrime);
				}
			}

			return primesValues;
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

			size_t possiblePrime = this->template getStartSieveNumber();
			WheelType _wheelSieve = this->template getWheelSieve();

			if (possiblePrime == this->template getBlockBeginNumber()) {
				primesFound = _wheelSieve.getNumberPrimesSievedByTheWheel();
			} else {
				primesFound = 0;
			}

			if (!(_wheelSieve.isNumberPossiblePrime(possiblePrime))) {
				possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime);
			}
			size_t maxRange = this->template getMaxRange();
			for (; possiblePrime <= maxRange; possiblePrime = _wheelSieve.getNextPossiblePrime(possiblePrime)) {
				if (!(this->template getPrimesBitsetValueMPI(possiblePrime))) {
					++primesFound;
				}
			}

			this->template setPrimesCount(primesFound);
			return primesFound;
		}

		// ------------------------------ <code from PrimesSieveParallelMultiplesOptimizedOpenMPI to avoid using virtuals in memory access functions (hotspots)> -------------------------------------------

};

