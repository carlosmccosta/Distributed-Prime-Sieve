#pragma once

#include "../PrimesSieve.h"

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

	public:
		PrimesSieveParallelMultiplesOptimizedOpenMP(size_t blockSizeInElements = 128 * 1024, size_t numberOfThreads = 0) :
				_blockSizeInElements(blockSizeInElements), _numberOfThreads(numberOfThreads) {
		}

		virtual ~PrimesSieveParallelMultiplesOptimizedOpenMP() {
		}

		void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template clearPrimesValues();
			this->template initPrimesBitSetSize(maxRange);

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			vector<pair<size_t, size_t> > sievingMultiples;

			this->template computeSievingPrimes(maxRangeSquareRoot, sievingMultiples);
			this->template removeComposites(maxRangeSquareRoot, maxRange, sievingMultiples);

			this->template getPerformanceTimer().stop();
		}

		virtual void computeSievingPrimes(size_t maxRangeSquareRoot, vector<pair<size_t, size_t> >& sievingMultiples) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumberOpenMP(maxRangeSquareRoot);
			size_t blockBeginNumber = getBlockBeginNumber();

			if (maxRangeSquareRoot < blockBeginNumber) {
				return;
			}

			size_t blockIndexBegin = this->template getBitsetPositionToNumberOpenMP(blockBeginNumber);
			size_t blockIndexEnd = min(blockIndexBegin + _blockSizeInElements, maxIndexRangeSquareRoot);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot, sievingMultiples);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;
				if (blockIndexEnd > maxIndexRangeSquareRoot) {
					blockIndexEnd = maxIndexRangeSquareRoot + 1;
				}

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
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

			#pragma omp parallel for \
			default(shared) \
			firstprivate(maxIndexRange, blockIndexSquareRoot, blockSizeInElements, sievingMultiples, priviousBlockNumber) \
			schedule(guided, 64) \
			num_threads(numberThreadsToUse)
			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + blockIndexSquareRoot;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				if (blockIndexEnd > maxIndexRange) {
					blockIndexEnd = maxIndexRange + 1;
				}

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPositionOpenMP(blockIndexEnd);

				if (blockNumber == 0) {
					sievingMultiples = sievingMultiplesFirstBlock;
				} else if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
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

		size_t getNumberOfThreads() const
		{
			return _numberOfThreads;
		}

		void setNumberOfThreads(size_t numberOfThreads)
		{
			_numberOfThreads = numberOfThreads;
		}
	};

