#pragma once

#include "../PrimesSieve.h"

#include <cmath>
#include <algorithm>
#include <utility>
#include <omp.h>

using std::sqrt;
using std::min;
using std::pair;

template<typename FlagsContainer>
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

			// adjustment of maxRange in order to calculate the primes <= maxRange instead of < maxRange
			if (maxRange % 2 == 0) {
				++maxRange;
			} else {
				maxRange += 2;
			}

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);

			this->template computeSievingPrimes(maxRangeSquareRoot);
			this->template removeComposites(maxRangeSquareRoot, maxRange);

			this->template getPerformanceTimer().stop();
		}

		virtual void computeSievingPrimes(size_t maxRangeSquareRoot) {
			size_t maxIndexRangeSquareRoot = this->template getBitsetPositionToNumber(maxRangeSquareRoot);
			size_t blockBeginNumber = getBlockBeginNumber();
			size_t blockIndexBegin = this->template getBitsetPositionToNumber(blockBeginNumber);
			size_t blockIndexEnd = min(_blockSizeInElements, maxIndexRangeSquareRoot);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

			if (blockEndNumber < blockBeginNumber) {
				return;
			}

			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);

			size_t numberBlocks = ceil((double) (maxRangeSquareRoot - blockBeginNumber) / (double) _blockSizeInElements);
			vector<pair<size_t, size_t> > sievingMultiples;

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;
				if (blockIndexEnd > maxIndexRangeSquareRoot) {
					blockIndexEnd = maxIndexRangeSquareRoot;
				}

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);
			}
		}

		virtual void removeComposites(size_t maxRangeSquareRoot, size_t maxRange) {
			const size_t blockSizeInElements = _blockSizeInElements;
			const size_t maxIndexRange = this->template getNumberBitsToStore(maxRange) - 1;
			const size_t blockIndexSquareRoot = this->template getBitsetPositionToNumber(maxRangeSquareRoot);

			const size_t numberBlocks = ceil((double) (maxIndexRange - blockIndexSquareRoot) / (double) blockSizeInElements);

			if (_numberOfThreads != 0) {
				if (numberBlocks < _numberOfThreads) {
					omp_set_num_threads(numberBlocks);
				} else {
					omp_set_num_threads(_numberOfThreads);
				}
			}

			vector<pair<size_t, size_t> > sievingMultiples;
			size_t priviousBlockNumber = -1;

//			#pragma omp threadprivate(sievingMultiples)
			#pragma omp parallel for \
			default(shared) \
			firstprivate(maxIndexRange, blockIndexSquareRoot, blockSizeInElements, sievingMultiples, priviousBlockNumber) \
			schedule(guided, 64)
			for (size_t blockNumber = 0; blockNumber < numberBlocks; ++blockNumber) {
				size_t blockIndexBegin = blockNumber * blockSizeInElements + blockIndexSquareRoot;
				size_t blockIndexEnd = blockIndexBegin + blockSizeInElements;

				if (blockIndexEnd > maxIndexRange) {
					blockIndexEnd = maxIndexRange;
				}

				size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				if (sievingMultiples.empty() || blockNumber != priviousBlockNumber + 1) {
					computeSievingMultiples(blockBeginNumber, blockEndNumber, sievingMultiples);
				}
				priviousBlockNumber = blockNumber;
				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, sievingMultiples);
			}
		}

		virtual void computeSievingMultiples(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, vector<pair<size_t, size_t> >& sievingMultiples) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot) = 0;
		virtual void initPrimesBitSetSize(size_t maxRange) = 0;

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
		}

		virtual size_t getBlockBeginNumber() {
			return this->template getNumberAssociatedWithBitsetPosition(0);
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

