#pragma once

#include "../PrimesSieve.h"

#include <cmath>
#include <algorithm>

using std::sqrt;
using std::min;

template<typename FlagsContainer>
class PrimesSieveSequencialMultiplesOptimized: public PrimesSieve<FlagsContainer> {
	protected:
		size_t _blockSizeInElements;

	public:
		PrimesSieveSequencialMultiplesOptimized(size_t blockSizeInElements = 128 * 1024) :
				_blockSizeInElements(blockSizeInElements) {
		}

		virtual ~PrimesSieveSequencialMultiplesOptimized() {
		}

		void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template clearPrimesValues();
			this->template initPrimesBitSetBlock(maxRange);

			// adjustment of maxRange in order to calculate the primes <= maxRange instead of < maxRange
			if (maxRange % 2 == 0) {
				++maxRange;
			} else {
				maxRange += 2;
			}

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			size_t maxIndexRange = this->template getNumberBitsToStore(maxRange) - 1;
			size_t numberBlocks = ceil((double) maxIndexRange / (double) _blockSizeInElements);
			size_t blockIndexBegin = 0;
			size_t blockIndexEnd = min(_blockSizeInElements, maxIndexRange);

			size_t blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);
			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;
				if (blockIndexEnd > maxIndexRange) {
					blockIndexEnd = maxIndexRange;
				}

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, blockIndexBegin);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);
			}

			this->template getPerformanceTimer().stop();
		}

		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot) = 0;
		virtual void initPrimesBitSetBlock(size_t maxRange) = 0;

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
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
	};

