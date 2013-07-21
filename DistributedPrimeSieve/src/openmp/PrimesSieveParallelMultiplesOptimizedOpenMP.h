#pragma once

#include "../PrimesSieve.h"
#include "../lib/PrimesUtils.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <omp.h>

using std::sqrt;
using std::min;
using std::pair;

template<typename FlagsContainer, typename WheelType>
class PrimesSieveParallelMultiplesOptimizedOpenMP: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;
		size_t _numberOfThreads;

		vector<size_t> _sievingPrimes;

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMP(size_t blockSizeInElements = 128 * 1024, size_t numberOfThreads = 0) :
				_blockSizeInElements(blockSizeInElements), _numberOfThreads(numberOfThreads) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMP() {
		}

		virtual void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template clearPrimesValues();

			PerformanceTimer memoryAllocationTimer;
			memoryAllocationTimer.start();
			this->template initPrimesBitSetSize(maxRange);
			memoryAllocationTimer.stop();
			size_t numberBitsToStore = this->template getNumberBitsToStore(maxRange);
			cout << "    --> Allocated memory for " << numberBitsToStore << " numbers in " << memoryAllocationTimer.getElapsedTimeFormated() << endl;

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			vector<pair<size_t, size_t> > sievingMultiples;

			// compute sieving primes
			PerformanceTimer computeSievingPrimesTimer;
			computeSievingPrimesTimer.start();
			this->template computeSievingPrimes(maxRangeSquareRoot, maxRange, sievingMultiples);
			computeSievingPrimesTimer.stop();
			cout << "    --> Computed " << sievingMultiples.size() << " sieving primes in " << computeSievingPrimesTimer.getElapsedTimeFormated() << endl;

			this->template removeComposites(maxRangeSquareRoot, maxRange, sievingMultiples);

			this->template getPerformanceTimer().stop();
		}

		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumberOpenMP(maxRangeSquareRoot);
			size_t blockBeginNumber = getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t blockIndexBegin = this->template getBitsetPositionToNumberOpenMP(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + _blockSizeInElements, maxIndexRangeSquareRoot + 1);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

				if (blockEndNumber >= maxRange) {
					blockEndNumber = maxRange + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocksParallel(blockBeginNumber, blockEndNumber, sievingMultiples);
//				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);
			}
		}

		virtual void removeComposites(size_t maxRangeSquareRoot, size_t maxRange, vector<pair<size_t, size_t> >& sievingMultiplesFirstBlock) {
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t maxIndexRange = this->template getNumberBitsToStore(maxRange) - 1;
			const size_t blockIndexSquareRoot = this->template getBitsetPositionToNumberOpenMP(maxRangeSquareRoot);

			const size_t numberBlocks = ceil((double) (maxIndexRange - blockIndexSquareRoot) / (double) blockSizeInElements);
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

#			pragma omp parallel for \
				if (numberBlocks > 32) \
				default(shared) \
				shared(sievingMultiplesFirstBlock) \
				firstprivate(maxIndexRange, blockIndexSquareRoot, blockSizeInElements, priviousBlockNumber, sievingMultiples) \
				schedule(guided, 64) \
				num_threads(numberThreadsToUse)
			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + blockIndexSquareRoot;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

				if (blockEndNumber > maxRange) {
					blockEndNumber = maxRange + 1;
				}

				if (blockNumber == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) {
			sievingMultiples.clear();
			size_t sievingPrimesSize = _sievingPrimes.size();
			for (size_t sievingPrimesIndex = 0; sievingPrimesIndex < sievingPrimesSize; ++sievingPrimesIndex) {
				size_t primeNumber = _sievingPrimes[sievingPrimesIndex];
				size_t primeMultiple = PrimesUtils::closestPrimeMultiple(primeNumber, blockBeginNumber);
				size_t primeMultipleIncrement = primeNumber << 1;

				if (primeMultiple < blockBeginNumber || primeMultiple == primeNumber) {
					primeMultiple += primeNumber;
				}

				if (primeMultiple % 2 == 0) {
					primeMultiple += primeNumber;
				}

				sievingMultiples.push_back(pair<size_t, size_t>(primeMultiple, primeMultipleIncrement));
			}
//					cout << "init sievingMultiples in block [" << blockBeginNumber << ", " << blockEndNumber << "]" << endl;
		}

		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeMultiplesOfPrimesFromPreviousBlocksParallel(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void initPrimesBitSetSize(size_t maxRange) = 0;
		virtual WheelType& getWheelSieve() = 0;

		virtual inline size_t getBitsetPositionToNumberOpenMP(size_t number) {
			return number;
		}

		virtual inline size_t getNumberAssociatedWithBitsetPositionOpenMP(size_t position) {
			return position;
		}

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
		}

		virtual size_t getBlockBeginNumber() {
			return this->template getNumberAssociatedWithBitsetPositionOpenMP(0);
		}

		inline size_t getBlockSizeInBytes() {
			return _blockSizeInElements / 8;
		}

		inline void setBlockSizeInBytes(size_t blockSize) {
			_blockSizeInElements = blockSize * 8;
		}

		inline size_t getBlockSizeInElements() const {
			return _blockSizeInElements;
		}

		inline void setBlockSizeInElements(size_t blockSizeInElements) {
			_blockSizeInElements = blockSizeInElements;
		}

		inline size_t getNumberOfThreads() const {
			return _numberOfThreads;
		}

		inline void setNumberOfThreads(size_t numberOfThreads) {
			_numberOfThreads = numberOfThreads;
		}

		inline vector<size_t>& getSievingPrimes() {
			return _sievingPrimes;
		}

		inline void setSievingPrimes(const vector<size_t>& sievingPrimes) {
			_sievingPrimes = sievingPrimes;
		}
	};

