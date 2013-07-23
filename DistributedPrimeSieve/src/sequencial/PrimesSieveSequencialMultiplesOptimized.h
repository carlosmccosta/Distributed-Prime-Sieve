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

		virtual void computePrimes(size_t maxRange) {
			this->template getPerformanceTimer().reset();
			this->template getPerformanceTimer().start();

			this->template clearPrimesValues();
			this->template initPrimesBitSetSize(maxRange);

			size_t maxRangeSquareRoot = (size_t) sqrt(maxRange);
			size_t maxIndexRange = this->template getNumberBitsToStore(maxRange);
			size_t numberBlocks = (size_t) ceil((double) maxIndexRange / (double) _blockSizeInElements);
			size_t blockIndexBegin = 0;
			size_t blockIndexEnd = min(_blockSizeInElements, maxIndexRange);

			size_t blockBeginNumber = getBlockBeginNumber();
			size_t blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);
			this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);

			for (size_t blockNumber = 1; blockNumber < numberBlocks; ++blockNumber) {
				blockIndexBegin = blockIndexEnd;
				blockIndexEnd += _blockSizeInElements;

				blockBeginNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexBegin);
				blockEndNumber = this->template getNumberAssociatedWithBitsetPosition(blockIndexEnd);

				if (blockEndNumber >= maxRange) {
					blockEndNumber = maxRange + 1;
				}

				this->template removeMultiplesOfPrimesFromPreviousBlocks(blockBeginNumber, blockEndNumber, blockIndexBegin);
				this->template calculatePrimesInBlock(blockBeginNumber, blockEndNumber, maxRangeSquareRoot);
			}

			this->template getPerformanceTimer().stop();
		}

		virtual void removeMultiplesOfPrimesFromPreviousBlocks(size_t blockBeginNumber, size_t blockEndNumber, size_t blockIndexBegin) = 0;
		virtual void calculatePrimesInBlock(size_t primeNumber, size_t maxNumberInBlock, size_t maxRangeSquareRoot) = 0;
		virtual void initPrimesBitSetSize(size_t maxRange) = 0;

		virtual void clearPrimesValues() {
			this->template getPrimesValues().clear();
		}

		virtual size_t getBlockBeginNumber() {
			return this->template getNumberAssociatedWithBitsetPosition(0);
		}

		virtual inline size_t getNumberAssociatedWithBitsetPosition(size_t position) {
			return this->PrimesSieve<FlagsContainer>::template getNumberAssociatedWithBitsetPosition(position);
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

